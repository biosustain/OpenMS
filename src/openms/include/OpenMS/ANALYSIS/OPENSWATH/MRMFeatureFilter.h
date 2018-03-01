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

#ifndef OPENMS_ANALYSIS_OPENSWATH_MRMFEATUREFILTER_H
#define OPENMS_ANALYSIS_OPENSWATH_MRMFEATUREFILTER_H

#include <OpenMS/ANALYSIS/MRM/ReactionMonitoringTransition.h>
#include <OpenMS/ANALYSIS/OPENSWATH/MRMFeatureQC.h>
#include <OpenMS/ANALYSIS/TARGETED/TargetedExperiment.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/FORMAT/QcMLFile.h>
#include <OpenMS/FORMAT/TraMLFile.h>
#include <OpenMS/KERNEL/Feature.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/KERNEL/MRMFeature.h>

namespace OpenMS
{
  /**
    @brief The MRMFeatureFilter either flags or filters out components
    and/or transitions that do not pass the QC criteria.

    @htmlinclude OpenMS_MRMFeatureFilter.parameters
  */
  class OPENMS_DLLAPI MRMFeatureFilter :
    public DefaultParamHandler
  {
public:
    /// Constructor
    MRMFeatureFilter();

    /// Destructor
    ~MRMFeatureFilter() = default;

    /// Synchronize members with param class
    void updateMembers_() override;

    enum C_Metrics : Size {RT = 0, INTENSITY, QUALITY, METAVALUE};
    enum CG_Metrics : Size {HEAVY = 0, LIGHT, DETECTING, QUANTIFYING, IDENTIFYING, TRANSITIONS, ION};

    /**
      @brief Flags or filters features and subordinates in a FeatureMap

      @param features FeatureMap to flag or filter
      @param filter_criteria MRMFeatureQC class defining QC parameters
      @param transitions transitions from a TargetedExperiment
    */
    void FilterFeatureMap(
      FeatureMap& features,
      const MRMFeatureQC& filter_criteria,
      const TargetedExperiment& transitions
    ) const;

    /*
      @brief Converts a FeatureMap to a qcMLFile::Attachment

      @param features FeatureMap to flag or filter
      @param attachment qcML Attachment

    */
    // void FeatureMapToAttachment(FeatureMap& features, QcMLFile::Attachment& attachment);

    /**
      @brief Calculates the ion ratio between two transitions

      @param component_1 component of the numerator
      @param component_2 component of the denominator
      @param feature_name name of the feature to calculate the ratio on
       e.g., peak_apex, peak_area

      @return The ratio.
    */
    double calculateIonRatio(const Feature& component_1, const Feature& component_2, const String& feature_name) const;

    /*
      @brief Calculates the retention time difference between two features

      @param component_1 First eluting component
      @param component_2 Second eluting component

      @return The difference.
    */
    // double calculateRTDifference(Feature& component_1, Feature& component_2);

    /*
      @brief Calculates the resolution between two features

      @param component_1 component 1
      @param component_2 component 2

      @return The difference.
    */
    // double calculateResolution(Feature& component_1, Feature& component_2);

    /**
      @brief Checks if the metaValue is within the user specified range

      @param component component of the numerator
      @param meta_value_key Name of the metaValue
      @param meta_value_l Lower bound (inclusive) for the metaValue range
      @param meta_value_u Upper bound (inclusive) for the metaValue range

      @return True if the metaValue is within the bounds, and False otherwise.
    */
    bool checkMetaValue(
      const Feature& component,
      const String& meta_value_key,
      const double& meta_value_l,
      const double& meta_value_u
    ) const;

    /**
      @brief Count the number of heavy/light labels and quantifying/detecting/identifying transitions

      @param component component_group with subordinates
      @param transitions transitions from a TargetedExperiment

      @return Map of labels/transition types and their corresponding number.
    */
    std::map<String, Int> countLabelsAndTransitionTypes(
      const Feature& component_group,
      const TargetedExperiment& transitions
    ) const;

    /**
      @brief Sorts, removes duplicates, and concatenates a list of Strings

      @param str_vec vector of Strings
      @param delim token to seperate Strings in the list

      @return A concatenated string.
    */
    String uniqueJoin(std::vector<String>& str_vec, const String& delim) const;

    /**
      @brief Tests whether the given value is in the range [value_l, value_u].

      @tparam T Type of the given inputs (e.g., `double`)

      @param[in] value Tested value
      @param[in] value_l Lower bound
      @param[in] value_u Upper bound

      @return `true` if the value is in the range `[value_l, value_u]`,
      `false` otherwise.
    */
    template <typename T>
    bool isValueWithinRange(const T& value, const T& value_l, const T& value_u) const;

protected:
    /// The number of QC metrics for a single component
    Size N_C_QC_METRICS = 4;
    /// The number of QC metrics for a component group
    Size N_CG_QC_METRICS = 7;
    // TODO document this
    bool allQCMetricsPass(std::vector<bool> qc_metrics_pass) const;
    // TODO document this
    double computeScore(std::vector<bool> qc_metrics_pass) const;

private:
    /// Flag or filter (i.e., remove) features that do not pass the QC
    String flag_or_filter_;
    /// include the data points for the extracted ion chromatogram (XIC) in the attachment
    bool report_xic_;
    /// include the data points for the total ion chromatogram (TIC) in the attachment
    bool report_tic_;
    // qcMLFile Attachment
    // QcMLFile::Attachment attachment_;
    // FeatureMap
    // FeatureMap features_;
    // component group/peptide/compound QCs
    // std::map<String,std::vector<QcMLFile::QualityParameter>> component_group_qc_report_;
    // component/transition QCs
    // std::map<String,std::vector<QcMLFile::QualityParameter>> component_qc_report_;
    // multi transition QCs
    // std::map<std::vector<String>,std::vector<QcMLFile::QualityParameter>> multi_component_group_qc_report_;
  };
}

#endif // OPENMS_ANALYSIS_OPENSWATH_MRMFEATUREFILTER_H
