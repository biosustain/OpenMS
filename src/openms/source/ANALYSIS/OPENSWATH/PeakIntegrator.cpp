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
    estimateBackground_(chromatogram, left, right);
  }

  void PeakIntegrator::estimateBackground(
    const MSSpectrum& spectrum,
    const double& left,
    const double& right
  )
  {
    estimateBackground_(spectrum, left, right);
  }

  void PeakIntegrator::integratePeak(
    const MSChromatogram& chromatogram,
    const double& left,
    const double& right
  )
  {
    integratePeak_(chromatogram, left, right);
  }

  void PeakIntegrator::integratePeak(
    const MSSpectrum& spectrum,
    const double& left,
    const double& right
  )
  {
    integratePeak_(spectrum, left, right);
  }

  double PeakIntegrator::simpson(MSChromatogram::ConstIterator it_begin, MSChromatogram::ConstIterator it_end) const
  {
    return simpson_(it_begin, it_end);
  }

  double PeakIntegrator::simpson(MSSpectrum::ConstIterator it_begin, MSSpectrum::ConstIterator it_end) const
  {
    return simpson_(it_begin, it_end);
  }

  void PeakIntegrator::calculatePeakShapeMetrics(
    const MSChromatogram& chromatogram,
    const double& left,
    const double& right,
    PeakShapeMetrics_& peakShapeMetrics
  )
  {
    calculatePeakShapeMetrics_(chromatogram, left, right, peakShapeMetrics);
  }

  void PeakIntegrator::calculatePeakShapeMetrics(
    const MSSpectrum& spectrum,
    const double& left,
    const double& right,
    PeakShapeMetrics_& peakShapeMetrics
  )
  {
    calculatePeakShapeMetrics_(spectrum, left, right, peakShapeMetrics);
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
