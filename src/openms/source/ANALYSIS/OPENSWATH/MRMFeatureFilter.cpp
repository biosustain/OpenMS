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
// $Authors: Douglas McCloskeyt, Pasquale Domenico Colaianni $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/OPENSWATH/MRMFeatureFilter.h>

namespace OpenMS
{
  MRMFeatureFilter::MRMFeatureFilter() :
    DefaultParamHandler("MRMFeatureFilter")
  {
    defaults_.setValue("flag_or_filter", "flag", "Flag or Filter (i.e., remove) Components or transitions that do not pass the QC.", ListUtils::create<String>("advanced"));
    defaults_.setValidStrings("flag_or_filter", ListUtils::create<String>("flag,filter"));

    //TODO:  future implementation for QcML reporting
    defaults_.setValue("report_xic", "false", "Embed an image of the XIC in the QC report.", ListUtils::create<String>("advanced"));
    defaults_.setValidStrings("report_xic", ListUtils::create<String>("true,false"));

    //TODO:  future implementation for QcML reporting
    defaults_.setValue("report_tic", "false", "Embed an image of the TIC in the QC report.", ListUtils::create<String>("advanced"));
    defaults_.setValidStrings("report_tic", ListUtils::create<String>("true,false"));

    // write defaults into Param object param_
    defaultsToParam_();
    updateMembers_();
  }

  void MRMFeatureFilter::updateMembers_()
  {
    flag_or_filter_ = (String)param_.getValue("flag_or_filter");
    report_xic_ = param_.getValue("report_xic").toBool();
    report_tic_ = param_.getValue("report_tic").toBool();
  }

  void MRMFeatureFilter::FilterFeatureMap(
    FeatureMap& features,
    const MRMFeatureQC& filter_criteria,
    const TargetedExperiment& transitions
  ) const
  {
    // initialize QC variables
    FeatureMap features_filtered;
    const String fail_message_delim = ";";
    for (Feature& f : features)
    {
      const String component_group_name = (String)f.getMetaValue("PeptideRef");
      std::map<String, Int> n_lbl_trans_type = countLabelsAndTransitionTypes(f, transitions);
      // initialize the new feature and subordinates
      std::vector<Feature> subordinates_filtered;
      std::vector<bool> cg_qc_passes(N_CG_QC_METRICS, true);
      std::vector<String> cg_qc_fail_message_vec;
      // iterate through each component/sub-feature
      for (Feature& subordinate : f.getSubordinates())
      {
        const String component_name = (String)subordinate.getMetaValue("native_id");
        std::vector<String> c_qc_fail_message_vec;
        // iterate through multi-feature/multi-sub-feature QCs/filters
        // iterate through component_groups
        for (const MRMFeatureQC::ComponentGroupQCs& cg_cq : filter_criteria.component_group_qcs)
        {
          if (cg_cq.component_group_name == component_group_name)
          {
            // labels and transition counts QC
            if (!isWithinRange(n_lbl_trans_type["n_heavy"], cg_cq.n_heavy_l, cg_cq.n_heavy_u))
            {
              cg_qc_passes[CG_Metrics::HEAVY] = false;
              cg_qc_fail_message_vec.push_back("n_heavy");
            }
            if (!isWithinRange(n_lbl_trans_type["n_light"], cg_cq.n_light_l, cg_cq.n_light_u))
            {
              cg_qc_passes[CG_Metrics::LIGHT] = false;
              cg_qc_fail_message_vec.push_back("n_light");
            }
            if (!isWithinRange(n_lbl_trans_type["n_detecting"], cg_cq.n_detecting_l, cg_cq.n_detecting_u))
            {
              cg_qc_passes[CG_Metrics::DETECTING] = false;
              cg_qc_fail_message_vec.push_back("n_detecting");
            }
            if (!isWithinRange(n_lbl_trans_type["n_quantifying"], cg_cq.n_quantifying_l, cg_cq.n_quantifying_u))
            {
              cg_qc_passes[CG_Metrics::QUANTIFYING] = false;
              cg_qc_fail_message_vec.push_back("n_quantifying");
            }
            if (!isWithinRange(n_lbl_trans_type["n_identifying"], cg_cq.n_identifying_l, cg_cq.n_identifying_u))
            {
              cg_qc_passes[CG_Metrics::IDENTIFYING] = false;
              cg_qc_fail_message_vec.push_back("n_identifying");
            }
            if (!isWithinRange(n_lbl_trans_type["n_transitions"], cg_cq.n_transitions_l, cg_cq.n_transitions_u))
            {
              cg_qc_passes[CG_Metrics::TRANSITIONS] = false;
              cg_qc_fail_message_vec.push_back("n_transitions");
            }

            // ion ratio QC
            for (const Feature& s : f.getSubordinates())
            {
              const String component_name2 = (String)s.getMetaValue("native_id");
              // find the ion ratio pair
              if (
                !cg_cq.ion_ratio_pair_name_1.empty() &&
                !cg_cq.ion_ratio_pair_name_2.empty() &&
                cg_cq.ion_ratio_pair_name_1 == component_name &&
                cg_cq.ion_ratio_pair_name_2 == component_name2
              )
              {
                const double ion_ratio = calculateIonRatio(subordinate, s, cg_cq.ion_ratio_feature_name);
                if (!isWithinRange(ion_ratio, cg_cq.ion_ratio_l, cg_cq.ion_ratio_u))
                {
                  cg_qc_passes[CG_Metrics::ION] = false;
                  cg_qc_fail_message_vec.push_back("ion_ratio_pair[" + component_name + "/" + component_name2 + "]");
                }
              }
            }
          }
        }
        std::vector<bool> c_qc_passes(N_C_QC_METRICS, true);
        // iterate through feature/sub-feature QCs/filters
        for (const MRMFeatureQC::ComponentQCs& c_qc : filter_criteria.component_qcs)
        {
          if (c_qc.component_name == component_name)
          {
            if (!isWithinRange((double)subordinate.getRT(), c_qc.retention_time_l, c_qc.retention_time_u))
            {
              c_qc_passes[C_Metrics::RT] = false;
              c_qc_fail_message_vec.push_back("retention_time");
            }
            if (!isWithinRange((double)subordinate.getIntensity(), c_qc.intensity_l, c_qc.intensity_u))
            {
              c_qc_passes[C_Metrics::INTENSITY] = false;
              c_qc_fail_message_vec.push_back("intensity");
            }
            if (!isWithinRange((double)subordinate.getOverallQuality(), c_qc.overall_quality_l, c_qc.overall_quality_u))
            {
              c_qc_passes[C_Metrics::QUALITY] = false;
              c_qc_fail_message_vec.push_back("overall_quality");
            }
            for (const std::pair<String, std::pair<double,double>>& kv : c_qc.meta_value_qc)
            {
              if (!checkMetaValue(subordinate, kv.first, kv.second.first, kv.second.second))
              {
                c_qc_passes[C_Metrics::METAVALUE] = false;
                c_qc_fail_message_vec.push_back("metaValue[" + kv.first + "]");
              }
            }
          }
        }
        // Copy or Flag passing/failing subordinates
        if (flag_or_filter_ == "filter")
        {
          if (allQCMetricsPass(c_qc_passes))
          {
            subordinates_filtered.push_back(subordinate);
          }
        }
        else if (flag_or_filter_ == "flag")
        {
          if (allQCMetricsPass(c_qc_passes))
          {
            subordinate.setMetaValue("QC_transition_pass", true);
            subordinate.setMetaValue("QC_transition_message", "");
          }
          else
          {
            subordinate.setMetaValue("QC_transition_pass", false);
            const String c_qc_fail_message = uniqueJoin(c_qc_fail_message_vec, fail_message_delim);
            subordinate.setMetaValue("QC_transition_message", c_qc_fail_message);
          }
        }
        subordinate.setMetaValue("QC_transition_score", computeScore(c_qc_passes));
      }

      // Copy or Flag passing/failing Features
      if (flag_or_filter_ == "filter")
      {
        if (allQCMetricsPass(cg_qc_passes) && subordinates_filtered.size())
        {
          Feature feature_filtered(f);
          feature_filtered.setSubordinates(subordinates_filtered);
          features_filtered.push_back(feature_filtered);
        }
      }
      else if (flag_or_filter_ == "flag")
      {
        if (allQCMetricsPass(cg_qc_passes))
        {
          f.setMetaValue("QC_transition_group_pass", true);
          f.setMetaValue("QC_transition_group_message", "");
        }
        else
        {
          f.setMetaValue("QC_transition_group_pass", false);
          const String cg_qc_fail_message = uniqueJoin(cg_qc_fail_message_vec, fail_message_delim);
          f.setMetaValue("QC_transition_group_message", cg_qc_fail_message);
        }
      }
      f.setMetaValue("QC_transition_group_score", computeScore(cg_qc_passes));
    }

    // replace with the filtered featureMap
    if (flag_or_filter_ == "filter")
    {
      features = features_filtered;
    }
  }

  std::map<String,Int> MRMFeatureFilter::countLabelsAndTransitionTypes(
    const Feature& component_group,
    const TargetedExperiment& transitions
  ) const
  {
    Int n_heavy{0}, n_light{0}, n_quant{0}, n_detect{0}, n_ident{0}, n_trans{0};
    std::map<String, Int> output;
    for (const Feature& subordinate : component_group.getSubordinates())
    {
      ++n_trans;
      const String label_type = (String)subordinate.getMetaValue("LabelType");
      if (label_type == "Heavy") ++n_heavy;
      else if (label_type == "Light") ++n_light;
      std::vector<ReactionMonitoringTransition>::const_iterator transition = find_if(
        transitions.getTransitions().begin(),
        transitions.getTransitions().end(),
        [&subordinate](const ReactionMonitoringTransition& tr)
        {
          return subordinate.getMetaValue("native_id") == tr.getNativeID();
        }
      );
      if (transition->isQuantifyingTransition()) ++n_quant;
      if (transition->isIdentifyingTransition()) ++n_ident;
      if (transition->isDetectingTransition()) ++n_detect;
    }
    output["n_heavy"] = n_heavy;
    output["n_light"] = n_light;
    output["n_quantifying"] = n_quant;
    output["n_identifying"] = n_ident;
    output["n_detecting"] = n_detect;
    output["n_transitions"] = n_trans;
    return output;
  }

  double MRMFeatureFilter::calculateIonRatio(
    const Feature& component_1,
    const Feature& component_2,
    const String& feature_name
  ) const
  {
    double ratio{0};
    if (component_1.metaValueExists(feature_name) && component_2.metaValueExists(feature_name))
    {
      ratio = (double)component_1.getMetaValue(feature_name) / (double)component_2.getMetaValue(feature_name);
    }
    else if (component_1.metaValueExists(feature_name))
    {
      LOG_WARN << "Warning: no ion pair found for transition_id " << component_1.getMetaValue("native_id") << ".";
      ratio = component_1.getMetaValue(feature_name);
    }
    else
    {
      LOG_INFO << "Feature metaValue " << feature_name << " not found for transition_ids " << component_1.getMetaValue("native_id") << " and " << component_2.getMetaValue("native_id") << ".";
    }
    return ratio;
  }

  bool MRMFeatureFilter::checkMetaValue(
    const Feature& component,
    const String& meta_value_key,
    const double& meta_value_l,
    const double& meta_value_u
  ) const
  {
    bool check = true;
    if (component.metaValueExists(meta_value_key))
    {
      check = isWithinRange((double)component.getMetaValue(meta_value_key), meta_value_l, meta_value_u);
    }
    else
    {
      LOG_WARN << "Warning: no metaValue found for transition_id " << component.getMetaValue("native_id") << " for metaValue key " << meta_value_key << "." << std::endl;
    }
    return check;
  }

  String MRMFeatureFilter::uniqueJoin(std::vector<String>& str_vec, const String& delim) const
  {
    // remove duplicates
    std::sort(str_vec.begin(), str_vec.end());
    str_vec.erase(std::unique(str_vec.begin(), str_vec.end()), str_vec.end());
    // concatenate
    String str_cat = "";
    for (const String& str : str_vec)
    {
      str_cat += str + delim;
    }
    // remove trailing delim
    if (!str_cat.empty() && !delim.empty())
    {
      str_cat.erase(str_cat.length() - delim.length());
    }
    return str_cat;
  }

  template <typename T>
  bool MRMFeatureFilter::isWithinRange(const T& value, const T& value_l, const T& value_u) const
  {
    return value >= value_l && value <= value_u;
  }

  //TODO: Future addition to allow for generating a QcML attachment for QC reporting
  // void MRMFeatureFilter::FeatureMapToAttachment(FeatureMap& features, QcMLFile::Attachment& attachment)
  // {
  //   //TODO
  // }

  bool MRMFeatureFilter::allQCMetricsPass(std::vector<bool> qc_metrics_pass) const
  {
    return std::all_of(qc_metrics_pass.begin(), qc_metrics_pass.end(), [](bool pass)
    {
      return pass;
    });
  }

  double MRMFeatureFilter::computeScore(std::vector<bool> qc_metrics_pass) const
  {
    const double score = std::accumulate(qc_metrics_pass.begin(), qc_metrics_pass.end(), 0);
    return score / qc_metrics_pass.size();
  }
}
