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
// $Maintainer: Douglas McCloskey, Pasquale Domenico Colaianni $
// $Authors: Douglas McCloskey, Pasquale Domenico Colaianni $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/OPENSWATH/PeakIntegrator.h>

namespace OpenMS
{
  PeakIntegrator::PeakIntegrator() :
    DefaultParamHandler("PeakIntegrator")
  {
    getDefaultParameters(defaults_);
    defaultsToParam_(); // write defaults into Param object param_
  }

  PeakIntegrator::~PeakIntegrator() {}

  void PeakIntegrator::estimateBackground(
    const MSChromatogram& chromatogram,
    const double& left,
    const double& right
  )
  {
    const double int_l = chromatogram.RTBegin(left)->getIntensity();
    const double int_r = (chromatogram.RTEnd(right)-1)->getIntensity();
    const double delta_int = int_r - int_l;
    const double delta_rt = (chromatogram.RTEnd(right)-1)->getRT() - chromatogram.RTBegin(left)->getRT();
    const double rt_min = int_r <= int_l ? (chromatogram.RTEnd(right)-1)->getRT() : chromatogram.RTBegin(left)->getRT();
    const double delta_int_apex = std::fabs(delta_int) * std::fabs(rt_min - peak_apex_rt_) / delta_rt;
    background_height_ = std::min(int_r, int_l) + delta_int_apex;
    double background = 0.0;
    if (baseline_type_ == "base_to_base")
    {
      if (integration_type_ == "trapezoid" || integration_type_ == "simpson")
      {
        // formula for calculating the background using the trapezoidal rule
        // background = intensity_min*delta_rt + 0.5*delta_int*delta_rt;
        background = delta_rt * (std::min(int_r, int_l) + 0.5 * std::fabs(delta_int));
      }
      else if (integration_type_ == "intensity_sum")
      {
        // calculate the background using the formula
        // y = mx + b where x = retention time, m = slope, b = left intensity
        // sign of delta_int will determine line direction
        // background += delta_int / delta_rt * (it->getRT() - left) + int_l;
        UInt n_points = 0;
        for (auto it=chromatogram.RTBegin(left); it!=chromatogram.RTEnd(right); ++it, ++n_points)
        {
          background += it->getRT();
        }
        background = (background - n_points * chromatogram.RTBegin(left)->getRT()) * delta_int / delta_rt + n_points * int_l;
      }
    }
    else if (baseline_type_ == "vertical_division")
    {
      if (integration_type_ == "trapezoid" || integration_type_ == "simpson")
      {
        background = delta_rt * std::min(int_r, int_l);
      }
      else if (integration_type_ == "intensity_sum")
      {
        UInt n_points = 0;
        for (auto it=chromatogram.RTBegin(left); it!=chromatogram.RTEnd(right); ++it, ++n_points)
          ;
        background = std::min(int_r, int_l) * n_points;
      }
    }
    background_area_ = background;
  }

  void PeakIntegrator::integratePeak(
    const MSChromatogram& chromatogram,
    const double& left,
    const double& right
  )
  {
    peak_area_ = 0.0;
    peak_height_ = -1.0;
    peak_apex_rt_ = -1.0;
    UInt n_points = 0;
    for (auto it=chromatogram.RTBegin(left); it!=chromatogram.RTEnd(right); ++it, ++n_points)
    {
      if (peak_height_ < it->getIntensity())
      {
        peak_height_ = it->getIntensity();
        peak_apex_rt_ = it->getRT();
      }
    }

    if (getIntegrationType() == "trapezoid")
    {
      for (auto it=chromatogram.RTBegin(left); it!=chromatogram.RTEnd(right)-1; ++it)
      {
        peak_area_ += ((it+1)->getRT() - it->getRT()) * ((it->getIntensity() + (it+1)->getIntensity()) / 2.0);
      }
    }
    else if (getIntegrationType() == "simpson")
    {
      if (n_points < 3)
      {
        LOG_DEBUG << std::endl << "Error in integratePeak: number of points must be >=3 for Simpson's rule" << std::endl;
        return;
      }
      if (n_points % 2)
      {
        peak_area_ = simpson(chromatogram.RTBegin(left), chromatogram.RTEnd(right));
      }
      else
      {
        double areas[4] = {};
        areas[0] = simpson(chromatogram.RTBegin(left), chromatogram.RTEnd(right) - 1);   // without last point
        areas[1] = simpson(chromatogram.RTBegin(left) + 1, chromatogram.RTEnd(right));   // without first point
        if (chromatogram.begin() <= chromatogram.RTBegin(left) - 1)
        {
          areas[2] = simpson(chromatogram.RTBegin(left) - 1, chromatogram.RTEnd(right)); // with one more point on the left
        }
        if (chromatogram.RTEnd(right) < chromatogram.end())
        {
          areas[3] = simpson(chromatogram.RTBegin(left), chromatogram.RTEnd(right) + 1); // with one more point on the right
        }
        UInt valids = 0;
        for (auto area : areas)
        {
          if (area)
          {
            peak_area_ += area;
            ++valids;
          }
        }
        peak_area_ /= valids;
      }
    }
    else
    {
      std::cout << std::endl << "WARNING: intensity_sum method is being used." << std::endl;
      for (auto it=chromatogram.RTBegin(left); it!=chromatogram.RTEnd(right); ++it)
      {
        peak_area_ += it->getIntensity();
      }
    }
  }

  double PeakIntegrator::simpson(
    MSChromatogram::ConstIterator it_begin,
    MSChromatogram::ConstIterator it_end
  ) const
  {
    double integral = 0.0;
    for (auto it=it_begin+1; it<it_end-1; it=it+2)
    {
      const double h = it->getRT() - (it-1)->getRT();
      const double k = (it+1)->getRT() - it->getRT();
      const double y_h = (it-1)->getIntensity();
      const double y_0 = it->getIntensity();
      const double y_k = (it+1)->getIntensity();
      integral += (1.0/6.0) * (h+k) * ((2.0-k/h)*y_h + (pow(h+k,2)/(h*k))*y_0 + (2.0-h/k)*y_k);
    }
    return integral;
  }

  PeakIntegrator::PeakShapeMetrics_ PeakIntegrator::calculatePeakShapeMetrics(
    const MSChromatogram& chromatogram,
    const double& left,
    const double& right
  ) const
  {
    // TODO remove the three following lines
    // peak_height_ = 965356;
    // peak_apex_rt_ = 2.7045;
    // background_height_ = 723.5;
    PeakShapeMetrics_ peakShapeMetrics;
    peakShapeMetrics.points_across_baseline = 0;

    const double start_intensity = std::max(chromatogram.RTBegin(left)->getIntensity() - background_height_, 0.0);
    const double end_intensity = std::max((chromatogram.RTEnd(right)-1)->getIntensity() - background_height_, 0.0);
    // TODO remove the following couts
    std::cout << std::endl << "int left: " << chromatogram.RTBegin(left)->getIntensity() << std::endl;
    std::cout << "int right: " << (chromatogram.RTEnd(right)-1)->getIntensity() << std::endl;
    std::cout << "bg_height: " << background_height_ << std::endl;
    std::cout << "start_intensity: " << start_intensity << std:: endl;
    std::cout << "end_intensity: " << end_intensity << std::endl;

    MSChromatogram::ConstIterator it_rtbegin_l = chromatogram.RTBegin(left);
    MSChromatogram::ConstIterator it_rtend_apex = chromatogram.RTEnd(peak_apex_rt_);
    MSChromatogram::ConstIterator it_rtend_r = chromatogram.RTEnd(right);

    peakShapeMetrics.start_time_at_5 = findRTAtPoint(it_rtbegin_l, it_rtend_apex, 0.05, true);
    peakShapeMetrics.start_time_at_10 = findRTAtPoint(it_rtbegin_l, it_rtend_apex, 0.1, true);
    peakShapeMetrics.start_time_at_50 = findRTAtPoint(it_rtbegin_l, it_rtend_apex, 0.5, true);
    peakShapeMetrics.end_time_at_5 = findRTAtPoint(it_rtend_apex, it_rtend_r, 0.05, false);
    peakShapeMetrics.end_time_at_10 = findRTAtPoint(it_rtend_apex, it_rtend_r, 0.1, false);
    peakShapeMetrics.end_time_at_50 = findRTAtPoint(it_rtend_apex, it_rtend_r, 0.5, false);

    for (auto it = it_rtbegin_l; it != it_rtend_r; ++it)
    {
      peakShapeMetrics.points_across_baseline++;
      if (std::max(it->getIntensity() - background_height_, 0.0) >= 0.5 * peak_height_)
      {
        peakShapeMetrics.points_across_half_height++;
      }
    }

    // peak widths
    peakShapeMetrics.width_at_5 = peakShapeMetrics.end_time_at_5 - peakShapeMetrics.start_time_at_5;
    peakShapeMetrics.width_at_10 = peakShapeMetrics.end_time_at_10 - peakShapeMetrics.start_time_at_10;
    peakShapeMetrics.width_at_50 = peakShapeMetrics.end_time_at_50 - peakShapeMetrics.start_time_at_50;
    peakShapeMetrics.total_width = right - left;
    peakShapeMetrics.slope_of_baseline = end_intensity - start_intensity;
    peakShapeMetrics.baseline_delta_2_height = peakShapeMetrics.slope_of_baseline / peak_height_;

    // other
    peakShapeMetrics.tailing_factor =
      peakShapeMetrics.width_at_5 /
      std::min(
        peak_apex_rt_ - peakShapeMetrics.start_time_at_5,
        peakShapeMetrics.end_time_at_5 - peak_apex_rt_
      );

    peakShapeMetrics.asymmetry_factor =
      std::min(
        peak_apex_rt_ - peakShapeMetrics.start_time_at_10,
        peakShapeMetrics.end_time_at_10 - peak_apex_rt_
      ) /
      std::max(
        peak_apex_rt_ - peakShapeMetrics.start_time_at_10,
        peakShapeMetrics.end_time_at_10 - peak_apex_rt_
      );

    return peakShapeMetrics;
  }

  double PeakIntegrator::findRTAtPoint(
    MSChromatogram::ConstIterator it_begin,
    MSChromatogram::ConstIterator it_end,
    const double perc,
    const bool is_start_time
  ) const
  {
    MSChromatogram::ConstIterator it = std::find_if(it_begin, it_end, [perc, this](ChromatogramPeak pt)
    {
      return pt.getIntensity() >= perc * peak_height_;
    });
    const double int_r = std::max(it->getIntensity() - background_height_, 0.0); //background-subtracted intensity
    const double int_l = std::max((it-1)->getIntensity() - background_height_, 0.0); //background-subtracted intensity of the previous point
    const double rt_r = it->getMZ();
    const double rt_l = (it-1)->getMZ();
    return rt_r - (int_r - int_l) * (rt_r - rt_l) / (is_start_time ? int_r - perc * peak_height_ : perc * peak_height_ - int_r);
  }

  double PeakIntegrator::getPeakArea() const
  {
    return peak_area_;
  }

  double PeakIntegrator::getPeakHeight() const
  {
    return peak_height_;
  }

  double PeakIntegrator::getPeakApexRT() const
  {
    return peak_apex_rt_;
  }

  double PeakIntegrator::getBackgroundHeight() const
  {
    return background_height_;
  }

  double PeakIntegrator::getBackgroundArea() const
  {
    return background_area_;
  }

  void PeakIntegrator::setIntegrationType(const String& integration_type)
  {
    integration_type_ = integration_type;
  }

  String PeakIntegrator::getIntegrationType() const
  {
    return integration_type_;
  }

  void PeakIntegrator::setBaselineType(const String& baseline_type)
  {
    baseline_type_ = baseline_type;
  }

  String PeakIntegrator::getBaselineType() const
  {
    return baseline_type_;
  }

  void PeakIntegrator::setPeakModel(const String& peak_model)
  {
    peak_model_ = peak_model;
  }

  String PeakIntegrator::getPeakModel() const
  {
    return peak_model_;
  }

  void PeakIntegrator::getDefaultParameters(Param& params)
  {
    params.clear();
    // TODO improve descriptions
    params.setValue("integration_type", "trapezoid", "Integration method to use.");
    params.setValidStrings("integration_type", ListUtils::create<String>("intensity_sum,simpson,trapezoid"));

    params.setValue("baseline_type", "vertical_division", "Type of baseline to use.");
    params.setValidStrings("baseline_type", ListUtils::create<String>("base_to_base,vertical_division"));

    params.setValue("peak_model", "none", "Peak model.");
    params.setValidStrings("peak_model", ListUtils::create<String>("none"));
  }

  void PeakIntegrator::updateMembers_()
  {
    integration_type_ = (String)param_.getValue("integration_type");
    baseline_type_ = (String)param_.getValue("baseline_type");
    peak_model_ = (String)param_.getValue("peak_model");
  }
}
