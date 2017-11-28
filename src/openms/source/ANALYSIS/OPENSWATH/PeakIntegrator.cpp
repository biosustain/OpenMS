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

  double PeakIntegrator::estimateBackground(
    const MSChromatogram& chromatogram,
    const double& left,
    const double& right
  )
  {
    const double intensity_l = chromatogram.RTBegin(left)->getIntensity();
    const double intensity_r = (chromatogram.RTEnd(right)-1)->getIntensity();
    const double delta_int = intensity_r - intensity_l;
    double background = 0.0;
    if (integration_type_ == "trapezoid" || integration_type_ == "simpson")
    {
      const double intensity_min = std::min(intensity_r, intensity_l);
      const double delta_int_abs = std::fabs(delta_int);
      for (auto it=chromatogram.RTBegin(left); it<chromatogram.RTEnd(right)-1; ++it)
      {
        const double delta_rt = (it+1)->getRT() - it->getRT();
        background += intensity_min * delta_rt + 0.5 * delta_int_abs * delta_rt;
      }
    }
    else
    {
      const double delta_input_rt = right - left;
      for (auto it=chromatogram.RTBegin(left); it!=chromatogram.RTEnd(right); ++it)
      {
        // calculate the background using the formula
        // y = mx + b where x = retention time, m = slope, b = left intensity
        // sign of delta_int will determine line direction
        background += delta_int / delta_input_rt * (it->getRT() - left) + intensity_l;
      }
    }
    return background;
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
