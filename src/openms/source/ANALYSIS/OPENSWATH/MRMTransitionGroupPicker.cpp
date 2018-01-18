// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/OPENSWATH/MRMTransitionGroupPicker.h>

#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/PeakPickerHiRes.h>
#include <OpenMS/FILTERING/NOISEESTIMATION/SignalToNoiseEstimatorMedian.h>
#include <OpenMS/FILTERING/SMOOTHING/SavitzkyGolayFilter.h>
#include <OpenMS/FILTERING/SMOOTHING/GaussFilter.h>

namespace OpenMS
{

  // Simple linear interpolation at point x between x0 and x1
  double lin_interpolate(double x, double x0, double x1, double y0, double y1)
  {
    double slope = (y1 - y0) / (x1 - x0);
    double delta_y = (x - x0) * slope;
    return y0 + delta_y;
  }

  MRMTransitionGroupPicker::MRMTransitionGroupPicker() :
    DefaultParamHandler("MRMTransitionGroupPicker")
  {
    defaults_.setValue("stop_after_feature", -1, "Stop finding after feature (ordered by intensity; -1 means do not stop).");
    defaults_.setValue("stop_after_intensity_ratio", 0.0001, "Stop after reaching intensity ratio");
    defaults_.setValue("min_peak_width", -1.0, "Minimal peak width (s), discard all peaks below this value (-1 means no action).", ListUtils::create<String>("advanced"));

    defaults_.setValue("peak_integration", "original", "Calculate the peak area and height either the smoothed or the raw chromatogram data.", ListUtils::create<String>("advanced"));
    defaults_.setValidStrings("peak_integration", ListUtils::create<String>("original,smoothed"));

    defaults_.setValue("background_subtraction", "none", "Try to apply a background subtraction to the peak (experimental). The background is estimated as the average noise at the peak boundaries (original) or at the exact left and right peak positions (exact).  The same original or smoothed chromatogram specified by peak_integration will be used for background estimation.", ListUtils::create<String>("advanced"));
    defaults_.setValidStrings("background_subtraction", ListUtils::create<String>("none,original,exact"));

    defaults_.setValue("recalculate_peaks", "false", "Tries to get better peak picking by looking at peak consistency of all picked peaks. Tries to use the consensus (median) peak border if theof variation within the picked peaks is too large.", ListUtils::create<String>("advanced"));
    defaults_.setValidStrings("recalculate_peaks", ListUtils::create<String>("true,false"));

    defaults_.setValue("use_precursors", "false", "Use precursor chromatogram for peak picking", ListUtils::create<String>("advanced"));
    defaults_.setValidStrings("use_precursors", ListUtils::create<String>("true,false"));

    defaults_.setValue("recalculate_peaks_max_z", 1.0, "Determines the maximal Z-Score (difference measured in standard deviations) that is considered too large for peak boundaries. If the Z-Score is above this value, the median is used for peak boundaries (default value 1.0).", ListUtils::create<String>("advanced"));

    defaults_.setValue("minimal_quality", -10000.0, "Only if compute_peak_quality is set, this parameter will not consider peaks below this quality threshold", ListUtils::create<String>("advanced"));

    defaults_.setValue("resample_boundary", 15.0, "For computing peak quality, how many extra seconds should be sample left and right of the actual peak", ListUtils::create<String>("advanced"));

    defaults_.setValue("compute_peak_quality", "false", "Tries to compute a quality value for each peakgroup and detect outlier transitions. The resulting score is centered around zero and values above 0 are generally good and below -1 or -2 are usually bad.", ListUtils::create<String>("advanced"));
    defaults_.setValidStrings("compute_peak_quality", ListUtils::create<String>("true,false"));
    
    defaults_.setValue("compute_peak_shape_metrics", "false", "Calulates various peak shape metrics (e.g., tailing) that can be used for downstream QC/QA.", ListUtils::create<String>("advanced"));
    defaults_.setValidStrings("compute_peak_shape_metrics", ListUtils::create<String>("true,false"));

    defaults_.insert("PeakPickerMRM:", PeakPickerMRM().getDefaults());

    // write defaults into Param object param_
    defaultsToParam_();
    updateMembers_();
  }

  MRMTransitionGroupPicker::~MRMTransitionGroupPicker()
  {
  }

  MRMTransitionGroupPicker& MRMTransitionGroupPicker::operator=(const MRMTransitionGroupPicker& rhs)
  {
    if (&rhs == this)
      return *this;

    // don't copy parameters

    return *this;
  }

  void MRMTransitionGroupPicker::updateMembers_()
  {
    stop_after_feature_ = (int)param_.getValue("stop_after_feature");
    stop_after_intensity_ratio_ = (double)param_.getValue("stop_after_intensity_ratio");
    peak_integration_ = param_.getValue("peak_integration");
    background_subtraction_ = param_.getValue("background_subtraction");
    recalculate_peaks_ = (bool)param_.getValue("recalculate_peaks").toBool();
    use_precursors_ = (bool)param_.getValue("use_precursors").toBool();
    recalculate_peaks_max_z_ = (double)param_.getValue("recalculate_peaks_max_z");
    compute_peak_quality_ = (bool)param_.getValue("compute_peak_quality").toBool();
    compute_peak_shape_metrics_ = (bool)param_.getValue("compute_peak_shape_metrics").toBool();
    min_qual_ = (double)param_.getValue("minimal_quality");
    min_peak_width_ = (double)param_.getValue("min_peak_width");
    resample_boundary_ = (double)param_.getValue("resample_boundary");

    picker_.setParameters(param_.copy("PeakPickerMRM:", true));
  }
  
  void MRMTransitionGroupPicker::calculateBgEstimationAverage_(const MSChromatogram& chromatogram,
      double best_left, double best_right, double & background, double & avg_noise_level)
  {
    // determine (in the chromatogram) the intensity at the left / right border
    MSChromatogram::const_iterator it = chromatogram.begin();
    int nr_points = 0;
    for (; it != chromatogram.end(); ++it)
    {
      if (it->getMZ() > best_left)
      {
        nr_points++;
        break;
      }
    }
    double intensity_left = it->getIntensity();
    for (; it != chromatogram.end(); ++it)
    {
      if (it->getMZ() > best_right)
      {
        break;
      }
      nr_points++;
    }
    if (it == chromatogram.begin() || nr_points < 1)
    {
      // something is fishy, the endpoint of the peak is the beginning of the chromatogram
      std::cerr << "Tried to calculate background but no points were found " << std::endl;
      return;
    }

    // decrease the iterator and the nr_points by one (because we went one too far)
    double intensity_right = (it--)->getIntensity();
    nr_points--;

    avg_noise_level = (intensity_right + intensity_left) / 2;
    background = avg_noise_level * nr_points;
  }

  void MRMTransitionGroupPicker::calculateBgEstimationExact_(const MSChromatogram& chromatogram,
      double best_left, double best_right, double peak_height, double & background, double & avg_noise_level)
  {
    double peak_apex_RT(0.0);
    double max_height = 0.0;
    for (MSChromatogram::ConstIterator it = chromatogram.RTBegin(best_left); it != chromatogram.RTEnd(best_right); ++it)
    {
      if (max_height < it->getIntensity())
      {
        max_height = it->getIntensity();
        peak_apex_RT = it->getRT();
      }
    }
    Param params;
    params.setValue("baseline_type", PeakIntegrator::BASELINE_TYPE_BASETOBASE);
    params.setValue("integration_type", PeakIntegrator::INTEGRATION_TYPE_INTENSITYSUM);
    PeakIntegrator pi;
    pi.setParameters(params);
    PeakIntegrator::PeakBackground pb = pi.estimateBackground(chromatogram, best_left, best_right, peak_apex_RT);
    background = pb.area;
    avg_noise_level = pb.height;
  }

  void MRMTransitionGroupPicker::findLargestPeak(std::vector<MSChromatogram >& picked_chroms, int& chr_idx, int& peak_idx)
  {
    double largest = 0.0;
    ChromatogramPeak largest_pos;
    for (Size k = 0; k < picked_chroms.size(); k++)
    {
      for (Size i = 0; i < picked_chroms[k].size(); i++)
      {
        if (picked_chroms[k][i].getIntensity() > largest)
        {
          largest = picked_chroms[k][i].getIntensity();
          chr_idx = (int)k;
          peak_idx = (int)i;
        }
      }
    }
  }
  
  void MRMTransitionGroupPicker::calculatePeakApexInt_(const MSChromatogram& chromatogram,
    double best_left, double best_right, 
    ConvexHull2D::PointArrayType & hull_points,
    double & intensity_sum, 
    double & intensity_integral,
    double & rt_sum,
    double & peak_apex_int,
    double & peak_apex_rt)
  {
    Param params;
    params.setValue("integration_type", PeakIntegrator::INTEGRATION_TYPE_INTENSITYSUM);
    PeakIntegrator pi;
    pi.setParameters(params);
    PeakIntegrator::PeakArea pa = pi.integratePeak(chromatogram, best_left, best_right);
    intensity_sum = pa.area;
    params.setValue("integration_type", PeakIntegrator::INTEGRATION_TYPE_TRAPEZOID);
    pi.setParameters(params);
    pa = pi.integratePeak(chromatogram, best_left, best_right);
    hull_points = pa.hull_points;
    intensity_integral = pa.area;
    peak_apex_int = pa.height;
    peak_apex_rt = pa.apex_pos;
    rt_sum = 0.0;
    for (MSChromatogram::ConstIterator it = chromatogram.RTBegin(best_left); it != chromatogram.RTEnd(best_right); ++it)
    {
      rt_sum += it->getRT();
    }
  }

  void MRMTransitionGroupPicker::calculatePeakShapeMetrics_(const MSChromatogram& chromatogram, 
    double best_left, double best_right, 
    double peak_height, double peak_apex_rt, double avg_noise_level,
    PeakShapeMetrics_ & peakShapeMetrics)
  {
    PeakIntegrator pi;
    PeakIntegrator::PeakShapeMetrics psm = pi.calculatePeakShapeMetrics(chromatogram, best_left, best_right, peak_height, peak_apex_rt);
    peakShapeMetrics.width_at_5 = psm.width_at_5;
    peakShapeMetrics.width_at_10 = psm.width_at_10;
    peakShapeMetrics.width_at_50 = psm.width_at_50;
    peakShapeMetrics.start_time_at_5 = psm.start_position_at_5;
    peakShapeMetrics.start_time_at_10 = psm.start_position_at_10;
    peakShapeMetrics.start_time_at_50 = psm.start_position_at_50;
    peakShapeMetrics.end_time_at_5 = psm.end_position_at_5;
    peakShapeMetrics.end_time_at_10 = psm.end_position_at_10;
    peakShapeMetrics.end_time_at_50 = psm.end_position_at_50;
    peakShapeMetrics.total_width = psm.total_width;
    peakShapeMetrics.tailing_factor = psm.tailing_factor;
    peakShapeMetrics.asymmetry_factor = psm.asymmetry_factor;
    peakShapeMetrics.baseline_delta_2_height = psm.baseline_delta_2_height;
    peakShapeMetrics.slope_of_baseline = psm.slope_of_baseline;
    peakShapeMetrics.points_across_baseline = psm.points_across_baseline;
    peakShapeMetrics.points_across_half_height = psm.points_across_half_height;
  }

}

