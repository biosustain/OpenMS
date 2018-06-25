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

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/ANALYSIS/OPENSWATH/TargetedSpectraExtractor.h>
#include <OpenMS/ANALYSIS/OPENSWATH/TransitionTSVFile.h>
#include <OpenMS/FORMAT/MSPGenericFile.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/FORMAT/TextFile.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

vector<pair<String,String>> expected_
{
  {"Pyruvate_174", "Propanoic acid, 2-(methoxyimino)-, trimethylsilyl ester"},
  {"Lactate_219", "Propanoic acid, 2-[(trimethylsilyl)oxy]-, trimethylsilyl ester"},
  {"Succinate_247", "Butanedioic acid, bis(trimethylsilyl) ester"},
  {"Fumarate_245", "2-Butenedioic acid (E)-, bis(trimethylsilyl) ester"},
  {"Malate_245", "Butanedioic acid, [(trimethylsilyl)oxy]-, bis(trimethylsilyl) ester"},
  {"Erythritol_307", "meso-Erythritol, tetrakis(trimethylsilyl) ether"},
  {"aKG_304", "Pentanedioic acid, 2-(methoxyimino)-, bis(trimethylsilyl) ester"},
  {"PEP_370", "2-Propenoic acid, 2-[(trimethylsilyl)oxy]-, anhydride with bis(trimethylsilyl) hydrogen phosphate"},
  {"GAP_328", "Phosphoric acid, 3-(methoxyimino)-2-[(trimethylsilyl)oxy]propyl bis(trimethylsilyl) ester, (.+/-.)-"},
  {"DHAP_400", "Phosphoric acid, 2-(methoxyimino)-3-[(trimethylsilyl)oxy]propyl bis(trimethylsilyl) ester"},
  {"G3P_445", "Phosphoric acid, bis(trimethylsilyl) 2,3-bis[(trimethylsilyl)oxy]propyl ester"},
  {"3PG_459", "3,5-Dioxa-4-phospha-2-silaoctan-8-oic acid, 2,2-dimethyl-4,7-bis[(trimethylsilyl)oxy]-, trimethylsilyl ester, 4-oxide"},
  {"Citrate_465", "1,2,3-Propanetricarboxylic acid, 2-[(trimethylsilyl)oxy]-, tris(trimethylsilyl) ester"},
  {"E4P_217", "Phosphorimidic acid, N-methoxy-, 4-oxo-2,3-bis[(trimethylsilyl)oxy]butyl bis(trimethylsilyl) ester, [R-(R*,R*)]-"},
  // {"ribulose-5-phosphate", "D-erythro-2-Pentulose, 1,3,4-tris-O-(trimethylsilyl)-, O-methyloxime, 5-[bis(trimethylsilyl) phosphate]"},
  {"R5P_160", "D-Ribose, 2,3,4-tris-O-(trimethylsilyl)-, O-methyloxime, 5-[bis(trimethylsilyl) phosphate]"},
  {"F6P_217", "d-Fructose, 1,3,4,5-tetrakis-O-(trimethylsilyl)-, o-methyloxime, 6-[bis(trimethylsilyl) phosphate]"},
  {"G6P_160", "d-Glucose, 2,3,4,5-tetrakis-O-(trimethylsilyl)-, o-methyloxime, 6-[bis(trimethylsilyl) phosphate]"},
  // {"6-phosphogluconate", "Gluconic acid, 2,3,4,5-tetrakis-O-(trimethylsilyl)-, trimethylsilyl ester, bis(trimethylsilyl) phosphate, D-"},
  {"S7P_357", "D-Altro-2-Heptulose, 1,3,4,5,6-pentakis-O-(trimethylsilyl)-, O-methyloxime, 7-[bis(trimethylsilyl) phosphate]"},
  {"Pyruvate", "Propanoic acid, 2-(methoxyimino)-, trimethylsilyl ester"},
  {"Lactate", "Propanoic acid, 2-[(trimethylsilyl)oxy]-, trimethylsilyl ester"},
  {"Succinate", "Butanedioic acid, bis(trimethylsilyl) ester"},
  {"Fumarate", "2-Butenedioic acid (E)-, bis(trimethylsilyl) ester"},
  {"Malate", "Butanedioic acid, [(trimethylsilyl)oxy]-, bis(trimethylsilyl) ester"},
  {"Erythritol", "meso-Erythritol, tetrakis(trimethylsilyl) ether"},
  {"aKG", "Pentanedioic acid, 2-(methoxyimino)-, bis(trimethylsilyl) ester"},
  {"PEP", "2-Propenoic acid, 2-[(trimethylsilyl)oxy]-, anhydride with bis(trimethylsilyl) hydrogen phosphate"},
  {"GAP", "Phosphoric acid, 3-(methoxyimino)-2-[(trimethylsilyl)oxy]propyl bis(trimethylsilyl) ester, (.+/-.)-"},
  {"DHAP", "Phosphoric acid, 2-(methoxyimino)-3-[(trimethylsilyl)oxy]propyl bis(trimethylsilyl) ester"},
  {"G3P", "Phosphoric acid, bis(trimethylsilyl) 2,3-bis[(trimethylsilyl)oxy]propyl ester"},
  {"3PG", "3,5-Dioxa-4-phospha-2-silaoctan-8-oic acid, 2,2-dimethyl-4,7-bis[(trimethylsilyl)oxy]-, trimethylsilyl ester, 4-oxide"},
  {"Citrate", "1,2,3-Propanetricarboxylic acid, 2-[(trimethylsilyl)oxy]-, tris(trimethylsilyl) ester"},
  {"E4P", "Phosphorimidic acid, N-methoxy-, 4-oxo-2,3-bis[(trimethylsilyl)oxy]butyl bis(trimethylsilyl) ester, [R-(R*,R*)]-"},
  {"R5P", "D-Ribose, 2,3,4-tris-O-(trimethylsilyl)-, O-methyloxime, 5-[bis(trimethylsilyl) phosphate]"},
  {"F6P", "d-Fructose, 1,3,4,5-tetrakis-O-(trimethylsilyl)-, o-methyloxime, 6-[bis(trimethylsilyl) phosphate]"},
  {"G6P", "d-Glucose, 2,3,4,5-tetrakis-O-(trimethylsilyl)-, o-methyloxime, 6-[bis(trimethylsilyl) phosphate]"},
  {"S7P", "D-Altro-2-Heptulose, 1,3,4,5,6-pentakis-O-(trimethylsilyl)-, O-methyloxime, 7-[bis(trimethylsilyl) phosphate]"},
};

bool matchIsValid(const String& spectrum_name, const String& match_name)
{
  vector<pair<String,String>>::const_iterator it = find_if(expected_.cbegin(), expected_.cend(),
    [&spectrum_name] (const pair<String,String>& p)
    {
      const String& expected_name = p.first;
      return spectrum_name == expected_name;
    });
  if (it == expected_.cend()) return false;
  return it->second == match_name;
}

void createTabularSpectra(const String& filename, const vector<MSSpectrum>& spectra)
{
  TextFile txt;
  Size highest_size { 0 };
  vector<double> highest_intensities(spectra.size(), 0);
  for (Size i = 0; i < spectra.size(); ++i)
  {
    if (spectra[i].size() > highest_size)
    {
      highest_size = spectra[i].size();
    }
    for (const Peak1D& p : spectra[i])
    {
      if (p.getIntensity() > highest_intensities[i])
      {
        highest_intensities[i] = p.getIntensity();
      }
    }
  }
  txt.addLine("unknown_mz\tunknown_int\texpected_mz\texpected_int\tbest_mz\tbest_int");
  for (Size i = 0; i < highest_size; ++i)
  {
    String line = "";
    for (Size j = 0; j < spectra.size(); ++j)
    {
      if (i < spectra[j].size())
      {
        line += to_string(spectra[j][i].getMZ()) + "\t" + to_string(spectra[j][i].getIntensity() / highest_intensities[j]) + "\t";
      }
      else
      {
        line += "\t\t";
      }
    }
    line.pop_back();
    txt.addLine(line);
  }
  txt.store(filename);
}

void outputFeatureMapToCsv(const String& filename, const FeatureMap& features)
{
  FeatureMap f = features;
  f.sortByPosition();
  TextFile txt;
  // txt.addLine("feature_idx\tretention_time\tprecursor_mz\tmatch_name\tmatch_score");
  txt.addLine("feature_idx\tretention_time\tmatch_name\tmatch_score");
  for (Size i = 0; i < f.size(); ++i)
  {
    String line = to_string(i) + "\t";

    line += to_string(f[i].getRT()) + "\t";

    // line += to_string(f[i].getMZ()) + "\t";

    String name = f[i].metaValueExists("spectral_library_name") ? f[i].getMetaValue("spectral_library_name") : "";
    line += name + "\t";

    double score = f[i].metaValueExists("spectral_library_score") ? static_cast<double>(f[i].getMetaValue("spectral_library_score")) : -1.0;
    line += to_string(score);

    txt.addLine(line);
  }
  txt.store(filename);
}

vector<MSSpectrum>::const_iterator findSpectrumByName(const vector<MSSpectrum>& spectra, const String& name)
{
  vector<MSSpectrum>::const_iterator it;
  it = std::find_if(spectra.cbegin(), spectra.cend(), [&name] (const MSSpectrum& s)
    {
      return s.getName() == name;
    });
  if (it == spectra.cend())
  {
    throw Exception::ElementNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, name);
  }
  return it;
}

START_TEST(TargetedSpectraExtractor, "$Id$")

/////////////////////////////////////////////////////////////
// Raw spectrum data acquired in DDA mode (i.e., product ion full spectrum scan)
// measured on a QTRAP 5500 corresponding to C-Aconitate
// taken from E. coli grown on glucose M9 during steady-state
// for flux analysis.

const vector<double> mz = {
  61.92, 68.88, 71.4, 79.56, 84.6, 84.72, 84.84, 84.96, 85.08, 85.2, 85.32,
  85.44, 85.68, 85.8, 85.92, 86.04, 86.16, 86.28, 86.4, 87.72, 87.96, 88.08,
  90.36, 94.44, 99.84, 100.8, 101.04, 101.88, 102, 102.96, 110.16, 110.88,
  111, 111.12, 111.24, 111.84, 111.96, 112.08, 112.2, 112.32, 112.44, 112.56,
  112.68, 114, 128.16, 128.4, 128.88, 129, 129.12, 129.84, 129.96, 130.08,
  130.2, 130.32, 130.44, 130.56, 132.12, 138, 139.08, 140.16, 144.12, 146.04,
  146.16, 156, 156.12, 156.36, 173.76, 174, 174.12, 174.24, 174.36, 174.6, 175.08
};
const vector<double> intensity = {
  6705.41660838088, 1676.35415209522, 1676.35415209522, 1676.35415209522, 3352.70830419044,
  5029.06245628566, 8381.7707604761, 53643.332867047, 51966.9787149518, 6705.41660838088,
  8381.7707604761, 1676.35415209522, 11734.4790646665, 25145.3122814283, 68730.520235904,
  112315.72819038, 6705.41660838088, 6705.41660838088, 3352.70830419044, 1676.35415209522,
  1676.35415209522, 1676.35415209522, 3352.70830419044, 1676.35415209522, 1676.35415209522,
  1676.35415209522, 5029.06245628566, 3352.70830419044, 3352.70830419044, 3352.70830419044,
  1676.35415209522, 5029.06245628566, 3352.70830419044, 5029.06245628566, 3352.70830419044,
  5029.06245628566, 18439.8956730474, 20116.2498251426, 5029.06245628566, 1676.35415209522,
  1676.35415209522, 3352.70830419044, 3352.70830419044, 3352.70830419044, 6705.41660838088,
  1676.35415209522, 3352.70830419044, 3352.70830419044, 6705.41660838088, 5029.06245628566,
  10058.1249125713, 31850.7288898092, 10058.1249125713, 1676.35415209522, 1676.35415209522,
  3352.70830419044, 1676.35415209522, 1676.35415209522, 1676.35415209522, 3352.70830419044,
  1676.35415209522, 3352.70830419044, 1676.35415209522, 1676.35415209522, 5029.06245628566,
  1676.35415209522, 1676.35415209522, 1676.35415209522, 6705.41660838088, 11734.4790646665,
  6705.41660838088, 1676.35415209522, 1676.35415209522
};
MSSpectrum s;
for (Size i = 0; i != mz.size(); ++i)
{
  s.push_back(Peak1D(mz[i],intensity[i]));
}
const MSSpectrum spectrum = s;

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

TargetedSpectraExtractor* ptr = nullptr;
TargetedSpectraExtractor* null_ptr = nullptr;
const String experiment_path = OPENMS_GET_TEST_DATA_PATH("TargetedSpectraExtractor_13C1_spectra0to100.mzML");
const String target_list_path = OPENMS_GET_TEST_DATA_PATH("TargetedSpectraExtractor_13CFlux_TraML.csv");
MzMLFile mzml;
MSExperiment experiment;
TransitionTSVFile tsv_reader;
TargetedExperiment targeted_exp;
mzml.load(experiment_path, experiment);
Param tsv_params = tsv_reader.getParameters();
tsv_params.setValue("retentionTimeInterpretation", "minutes");
tsv_reader.setParameters(tsv_params);
tsv_reader.convertTSVToTargetedExperiment(target_list_path.c_str(), FileTypes::CSV, targeted_exp);

START_SECTION(TargetedSpectraExtractor())
{
  ptr = new TargetedSpectraExtractor();
  TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION(~TargetedSpectraExtractor())
{
  delete ptr;
}
END_SECTION

START_SECTION(const Param& getParameters() const)
{
  TargetedSpectraExtractor tse;
  Param params = tse.getParameters();
  TEST_EQUAL(params.getValue("rt_window"), 30.0)
  TEST_EQUAL(params.getValue("min_select_score"), 0.7)
  TEST_EQUAL(params.getValue("mz_tolerance"), 0.1)
  TEST_EQUAL(params.getValue("mz_unit_is_Da"), "true")
  TEST_EQUAL(params.getValue("SavitzkyGolayFilter:frame_length"), 15)
  TEST_EQUAL(params.getValue("SavitzkyGolayFilter:polynomial_order"), 3)
  TEST_EQUAL(params.getValue("GaussFilter:gaussian_width"), 0.2)
  TEST_EQUAL(params.getValue("use_gauss"), "true")
  TEST_EQUAL(params.getValue("PeakPickerHiRes:signal_to_noise"), 1.0)
  TEST_EQUAL(params.getValue("peak_height_min"), 0.0)
  TEST_EQUAL(params.getValue("peak_height_max"), std::numeric_limits<double>::max())
  TEST_EQUAL(params.getValue("fwhm_threshold"), 0.0)
  TEST_EQUAL(params.getValue("tic_weight"), 1.0)
  TEST_EQUAL(params.getValue("fwhm_weight"), 1.0)
  TEST_EQUAL(params.getValue("snr_weight"), 1.0)
  TEST_EQUAL(params.getValue("similarity_function"), "BinnedSpectralContrastAngle")
  TEST_EQUAL(params.getValue("top_matches_to_report"), 5)
  // TEST_EQUAL(params.getValue("bin_size"), 1.0)
  // TEST_EQUAL(params.getValue("peak_spread"), 0.0)
  // TEST_EQUAL(params.getValue("bin_offset"), 0.4)
  TEST_EQUAL(params.getValue("min_match_score"), 0.8)
}
END_SECTION

START_SECTION(void getDefaultParameters(Param& params) const)
{
  TargetedSpectraExtractor tse;
  Param params;
  tse.getDefaultParameters(params);
  TEST_EQUAL(params.getValue("rt_window"), 30.0)
  TEST_EQUAL(params.getValue("min_select_score"), 0.7)
  TEST_EQUAL(params.getValue("mz_tolerance"), 0.1)
  TEST_EQUAL(params.getValue("mz_unit_is_Da"), "true")
  TEST_EQUAL(params.getValue("use_gauss"), "true")
  TEST_EQUAL(params.getValue("peak_height_min"), 0.0)
  TEST_EQUAL(params.getValue("peak_height_max"), std::numeric_limits<double>::max())
  TEST_EQUAL(params.getValue("fwhm_threshold"), 0.0)
  TEST_EQUAL(params.getValue("tic_weight"), 1.0)
  TEST_EQUAL(params.getValue("fwhm_weight"), 1.0)
  TEST_EQUAL(params.getValue("snr_weight"), 1.0)
  TEST_EQUAL(params.getValue("similarity_function"), "BinnedSpectralContrastAngle")
  TEST_EQUAL(params.getValue("top_matches_to_report"), 5)
  // TEST_EQUAL(params.getValue("bin_size"), 1.0)
  // TEST_EQUAL(params.getValue("peak_spread"), 0.0)
  // TEST_EQUAL(params.getValue("bin_offset"), 0.4)
  TEST_EQUAL(params.getValue("min_match_score"), 0.8)
}
END_SECTION

START_SECTION(void annotateSpectra(
  const std::vector<MSSpectrum>& spectra,
  const TargetedExperiment& targeted_exp,
  std::vector<MSSpectrum>& annotated_spectra,
  FeatureMap& features,
  const bool compute_features = true
) const)
{
  TargetedSpectraExtractor tse;
  Param params = tse.getParameters();
  params.setValue("GaussFilter:gaussian_width", 0.25);
  params.setValue("peak_height_min", 15000.0);
  params.setValue("peak_height_max", 110000.0);
  params.setValue("fwhm_threshold", 0.23);
  tse.setParameters(params);

  const vector<MSSpectrum>& spectra = experiment.getSpectra();
  vector<MSSpectrum> annotated_spectra;
  FeatureMap features;

  tse.annotateSpectra(spectra, targeted_exp, annotated_spectra, features);

  TEST_EQUAL(annotated_spectra.size(), 30)
  TEST_EQUAL(annotated_spectra.size(), features.size())

  TEST_EQUAL(annotated_spectra[0].getName(), "met-L.met-L_m0-0")
  TEST_EQUAL(annotated_spectra[0].size(), 121)
  TEST_EQUAL(annotated_spectra[4].getName(), "met-L.met-L_m1-0")
  TEST_EQUAL(annotated_spectra[4].size(), 135)
  TEST_EQUAL(annotated_spectra[8].getName(), "asp-L.asp-L_m0-0")
  TEST_EQUAL(annotated_spectra[8].size(), 55)
  TEST_EQUAL(annotated_spectra[12].getName(), "asp-L.asp-L_m1-0")
  TEST_EQUAL(annotated_spectra[12].size(), 389)
  TEST_EQUAL(annotated_spectra[16].getName(), "asp-L.asp-L_m2-1")
  TEST_EQUAL(annotated_spectra[16].size(), 143)
  TEST_EQUAL(annotated_spectra[20].getName(), "glu-L.glu-L_m5-5")
  TEST_EQUAL(annotated_spectra[20].size(), 82)
  TEST_EQUAL(annotated_spectra[24].getName(), "glu-L.glu-L_m2-2")
  TEST_EQUAL(annotated_spectra[24].size(), 94)
  TEST_EQUAL(annotated_spectra[29].getName(), "skm.skm_m4-4")
  TEST_EQUAL(annotated_spectra[29].size(), 552)

  TEST_EQUAL(features[0].getMetaValue("transition_name"), "met-L.met-L_m0-0")
  TEST_REAL_SIMILAR(features[0].getRT(), 80.22100000002)
  TEST_REAL_SIMILAR(features[0].getMZ(), 148.052001953125)
  TEST_EQUAL(features[4].getMetaValue("transition_name"), "met-L.met-L_m1-0")
  TEST_REAL_SIMILAR(features[4].getRT(), 87.927)
  TEST_REAL_SIMILAR(features[4].getMZ(), 149.054992675781)
  TEST_EQUAL(features[8].getMetaValue("transition_name"), "asp-L.asp-L_m0-0")
  TEST_REAL_SIMILAR(features[8].getRT(), 126.37699999998)
  TEST_REAL_SIMILAR(features[8].getMZ(), 132.029998779297)
  TEST_EQUAL(features[12].getMetaValue("transition_name"), "asp-L.asp-L_m1-0")
  TEST_REAL_SIMILAR(features[12].getRT(), 131.73100000002)
  TEST_REAL_SIMILAR(features[12].getMZ(), 133.033004760742)
  TEST_EQUAL(features[16].getMetaValue("transition_name"), "asp-L.asp-L_m2-1")
  TEST_REAL_SIMILAR(features[16].getRT(), 138.29599999998)
  TEST_REAL_SIMILAR(features[16].getMZ(), 134.035995483398)
  TEST_EQUAL(features[20].getMetaValue("transition_name"), "glu-L.glu-L_m5-5")
  TEST_REAL_SIMILAR(features[20].getRT(), 141.70399999998)
  TEST_REAL_SIMILAR(features[20].getMZ(), 151.061996459961)
  TEST_EQUAL(features[24].getMetaValue("transition_name"), "glu-L.glu-L_m2-2")
  TEST_REAL_SIMILAR(features[24].getRT(), 148.473)
  TEST_REAL_SIMILAR(features[24].getMZ(), 148.052001953125)
  TEST_EQUAL(features[29].getMetaValue("transition_name"), "skm.skm_m4-4")
  TEST_REAL_SIMILAR(features[29].getRT(), 166.95400000002)
  TEST_REAL_SIMILAR(features[29].getMZ(), 177.057998657227)
}
END_SECTION

START_SECTION(void annotateSpectra(
  const std::vector<MSSpectrum>& spectra,
  const TargetedExperiment& targeted_exp,
  std::vector<MSSpectrum>& annotated_spectra
) const)
{
  TargetedSpectraExtractor tse;
  Param params = tse.getParameters();
  params.setValue("GaussFilter:gaussian_width", 0.25);
  params.setValue("peak_height_min", 15000.0);
  params.setValue("peak_height_max", 110000.0);
  params.setValue("fwhm_threshold", 0.23);
  tse.setParameters(params);

  const vector<MSSpectrum>& spectra = experiment.getSpectra();
  vector<MSSpectrum> annotated_spectra;

  tse.annotateSpectra(spectra, targeted_exp, annotated_spectra);

  TEST_EQUAL(annotated_spectra.size(), 30)

  TEST_EQUAL(annotated_spectra[0].getName(), "met-L.met-L_m0-0")
  TEST_EQUAL(annotated_spectra[0].size(), 121)
  TEST_EQUAL(annotated_spectra[4].getName(), "met-L.met-L_m1-0")
  TEST_EQUAL(annotated_spectra[4].size(), 135)
  TEST_EQUAL(annotated_spectra[20].getName(), "glu-L.glu-L_m5-5")
  TEST_EQUAL(annotated_spectra[20].size(), 82)
  TEST_EQUAL(annotated_spectra[24].getName(), "glu-L.glu-L_m2-2")
  TEST_EQUAL(annotated_spectra[24].size(), 94)
  TEST_EQUAL(annotated_spectra[29].getName(), "skm.skm_m4-4")
  TEST_EQUAL(annotated_spectra[29].size(), 552)
}
END_SECTION

START_SECTION(void pickSpectrum(const MSSpectrum& spectrum, MSSpectrum& picked_spectrum) const)
{
  MSSpectrum picked_spectrum;
  TargetedSpectraExtractor tse;
  Param params = tse.getParameters();
  params.setValue("GaussFilter:gaussian_width", 0.25);
  params.setValue("peak_height_min", 0.0);
  params.setValue("peak_height_max", 200000.0);
  params.setValue("fwhm_threshold", 0.0);
  tse.setParameters(params);

  tse.pickSpectrum(spectrum, picked_spectrum);

  TEST_NOT_EQUAL(spectrum.size(), picked_spectrum.size())
  TEST_EQUAL(picked_spectrum.size(), 6)
  MSSpectrum::Iterator it = picked_spectrum.begin();
  TEST_REAL_SIMILAR(it->getMZ(), 85.014)
  TEST_REAL_SIMILAR(it->getIntensity(), 60754.7)
  ++it;
  TEST_REAL_SIMILAR(it->getMZ(), 86.0196)
  TEST_REAL_SIMILAR(it->getIntensity(), 116036.0)
  ++it;
  TEST_REAL_SIMILAR(it->getMZ(), 112.033)
  TEST_REAL_SIMILAR(it->getIntensity(), 21941.9)
  ++it;
  TEST_REAL_SIMILAR(it->getMZ(), 129.396)
  TEST_REAL_SIMILAR(it->getIntensity(), 10575.5)
  ++it;
  TEST_REAL_SIMILAR(it->getMZ(), 130.081)
  TEST_REAL_SIMILAR(it->getIntensity(), 31838.1)
  ++it;
  TEST_REAL_SIMILAR(it->getMZ(), 174.24)
  TEST_REAL_SIMILAR(it->getIntensity(), 11731.3)

  params.setValue("peak_height_min", 15000.0);
  params.setValue("peak_height_max", 110000.0);
  tse.setParameters(params);

  tse.pickSpectrum(spectrum, picked_spectrum);

  // With the new filters on peaks' heights, less peaks get picked.
  TEST_EQUAL(picked_spectrum.size(), 3)
  it = picked_spectrum.begin();
  TEST_REAL_SIMILAR(it->getMZ(), 85.014)
  TEST_REAL_SIMILAR(it->getIntensity(), 60754.7)
  ++it;
  TEST_REAL_SIMILAR(it->getMZ(), 112.033)
  TEST_REAL_SIMILAR(it->getIntensity(), 21941.9)
  ++it;
  TEST_REAL_SIMILAR(it->getMZ(), 130.081)
  TEST_REAL_SIMILAR(it->getIntensity(), 31838.1)

  params.setValue("fwhm_threshold", 0.23);
  tse.setParameters(params);

  tse.pickSpectrum(spectrum, picked_spectrum);

  // Filtering also on fwhm, even less peaks get picked.
  TEST_EQUAL(picked_spectrum.size(), 2)
  it = picked_spectrum.begin();
  TEST_REAL_SIMILAR(it->getMZ(), 85.014)
  TEST_REAL_SIMILAR(it->getIntensity(), 60754.7)
  ++it;
  TEST_REAL_SIMILAR(it->getMZ(), 112.033)
  TEST_REAL_SIMILAR(it->getIntensity(), 21941.9)

  MSSpectrum unordered;
  unordered.emplace_back(Peak1D(10.0, 100.0));
  unordered.emplace_back(Peak1D(9.0, 100.0));
  TEST_EXCEPTION(Exception::IllegalArgument, tse.pickSpectrum(unordered, picked_spectrum));
}
END_SECTION

START_SECTION(void scoreSpectra(
  const std::vector<MSSpectrum>& annotated_spectra,
  const std::vector<MSSpectrum>& picked_spectra,
  FeatureMap& features,
  std::vector<MSSpectrum>& scored_spectra,
  const bool compute_features = true
) const)
{
  TargetedSpectraExtractor tse;
  Param params = tse.getParameters();
  params.setValue("GaussFilter:gaussian_width", 0.25);
  params.setValue("peak_height_min", 15000.0);
  params.setValue("peak_height_max", 110000.0);
  params.setValue("fwhm_threshold", 0.23);
  tse.setParameters(params);

  vector<MSSpectrum> annotated_spectra;
  FeatureMap features;
  const vector<MSSpectrum>& spectra = experiment.getSpectra();

  tse.annotateSpectra(spectra, targeted_exp, annotated_spectra, features);

  vector<MSSpectrum> picked_spectra(annotated_spectra.size());
  for (Size i = 0; i < annotated_spectra.size(); ++i)
  {
    tse.pickSpectrum(annotated_spectra[i], picked_spectra[i]);
  }

  for (Int i = annotated_spectra.size() - 1; i >= 0; --i)
  {
    if (picked_spectra[i].empty())
    {
      annotated_spectra.erase(annotated_spectra.begin() + i);
      picked_spectra.erase(picked_spectra.begin() + i);
      features.erase(features.begin() + i);
    }
  }
  TEST_EQUAL(annotated_spectra.size(), 20)
  TEST_EQUAL(annotated_spectra.size(), features.size())
  TEST_EQUAL(picked_spectra.size(), features.size())

  vector<MSSpectrum> scored_spectra;
  tse.scoreSpectra(annotated_spectra, picked_spectra, features, scored_spectra);

  TEST_EQUAL(scored_spectra.size(), 20)
  TEST_EQUAL(scored_spectra.size(), annotated_spectra.size())
  TEST_EQUAL(scored_spectra.size(), features.size())

  TEST_EQUAL(scored_spectra[0].getName(), "met-L.met-L_m0-0")
  TEST_REAL_SIMILAR(scored_spectra[0].getFloatDataArrays()[1][0], 15.2046270370483) // score
  TEST_REAL_SIMILAR(scored_spectra[0].getFloatDataArrays()[2][0], 5.3508939743042)  // total tic
  TEST_REAL_SIMILAR(scored_spectra[0].getFloatDataArrays()[3][0], 3.96267318725586) // inverse average fwhm
  TEST_REAL_SIMILAR(scored_spectra[0].getFloatDataArrays()[4][0], 5.89106035232544) // average snr

  TEST_EQUAL(scored_spectra[4].getName(), "asp-L.asp-L_m1-0")
  TEST_REAL_SIMILAR(scored_spectra[4].getFloatDataArrays()[1][0], 10.8893)
  TEST_REAL_SIMILAR(scored_spectra[4].getFloatDataArrays()[2][0], 6.49946)
  TEST_REAL_SIMILAR(scored_spectra[4].getFloatDataArrays()[3][0], 2.65215)
  TEST_REAL_SIMILAR(scored_spectra[4].getFloatDataArrays()[4][0], 1.73772)

  TEST_EQUAL(scored_spectra[8].getName(), "asp-L.asp-L_m2-1")
  TEST_REAL_SIMILAR(scored_spectra[8].getFloatDataArrays()[1][0], 16.1929)
  TEST_REAL_SIMILAR(scored_spectra[8].getFloatDataArrays()[2][0], 5.52142)
  TEST_REAL_SIMILAR(scored_spectra[8].getFloatDataArrays()[3][0], 3.44492)
  TEST_REAL_SIMILAR(scored_spectra[8].getFloatDataArrays()[4][0], 7.22662)

  TEST_EQUAL(scored_spectra[11].getName(), "asp-L.asp-L_m2-2")
  TEST_REAL_SIMILAR(scored_spectra[11].getFloatDataArrays()[1][0], 17.4552)
  TEST_REAL_SIMILAR(scored_spectra[11].getFloatDataArrays()[2][0], 5.48532)
  TEST_REAL_SIMILAR(scored_spectra[11].getFloatDataArrays()[3][0], 3.78555)
  TEST_REAL_SIMILAR(scored_spectra[11].getFloatDataArrays()[4][0], 8.18436)

  TEST_EQUAL(scored_spectra[15].getName(), "glu-L.glu-L_m1-1")
  TEST_REAL_SIMILAR(scored_spectra[15].getFloatDataArrays()[1][0], 13.5799)
  TEST_REAL_SIMILAR(scored_spectra[15].getFloatDataArrays()[2][0], 5.49089)
  TEST_REAL_SIMILAR(scored_spectra[15].getFloatDataArrays()[3][0], 3.53584)
  TEST_REAL_SIMILAR(scored_spectra[15].getFloatDataArrays()[4][0], 4.55314)

  TEST_EQUAL(scored_spectra[19].getName(), "skm.skm_m4-4")
  TEST_REAL_SIMILAR(scored_spectra[19].getFloatDataArrays()[1][0], 10.5746)
  TEST_REAL_SIMILAR(scored_spectra[19].getFloatDataArrays()[2][0], 6.60354)
  TEST_REAL_SIMILAR(scored_spectra[19].getFloatDataArrays()[3][0], 2.02869)
  TEST_REAL_SIMILAR(scored_spectra[19].getFloatDataArrays()[4][0], 1.94236)

  TEST_EQUAL(features[0].getMetaValue("transition_name"), "met-L.met-L_m0-0")
  TEST_REAL_SIMILAR(features[0].getIntensity(), 15.2046270370483)                  // score
  TEST_REAL_SIMILAR(features[0].getMetaValue("log10_total_tic"), 5.3508939743042)  // total tic
  TEST_REAL_SIMILAR(features[0].getMetaValue("inverse_avgFWHM"), 3.96267318725586) // inverse average fwhm
  TEST_REAL_SIMILAR(features[0].getMetaValue("avgSNR"), 5.89106035232544)          // average snr
  TEST_REAL_SIMILAR(features[0].getMetaValue("avgFWHM"), 0.252354895075162)        // average fwhm

  TEST_EQUAL(features[4].getMetaValue("transition_name"), "asp-L.asp-L_m1-0")
  TEST_REAL_SIMILAR(features[4].getIntensity(), 10.8893)
  TEST_REAL_SIMILAR(features[4].getMetaValue("log10_total_tic"), 6.49945796336373)
  TEST_REAL_SIMILAR(features[4].getMetaValue("inverse_avgFWHM"), 2.65214624318674)
  TEST_REAL_SIMILAR(features[4].getMetaValue("avgSNR"), 1.73772000291411)
  TEST_REAL_SIMILAR(features[4].getMetaValue("avgFWHM"), 0.377053114084097)

  TEST_EQUAL(features[8].getMetaValue("transition_name"), "asp-L.asp-L_m2-1")
  TEST_REAL_SIMILAR(features[8].getIntensity(), 16.1929)
  TEST_REAL_SIMILAR(features[8].getMetaValue("log10_total_tic"), 5.52141560620828)
  TEST_REAL_SIMILAR(features[8].getMetaValue("inverse_avgFWHM"), 3.44491858720322)
  TEST_REAL_SIMILAR(features[8].getMetaValue("avgSNR"), 7.22661551261844)
  TEST_REAL_SIMILAR(features[8].getMetaValue("avgFWHM"), 0.290282621979713)

  TEST_EQUAL(features[11].getMetaValue("transition_name"), "asp-L.asp-L_m2-2")
  TEST_REAL_SIMILAR(features[11].getIntensity(), 17.4552)
  TEST_REAL_SIMILAR(features[11].getMetaValue("log10_total_tic"), 5.48531541983726)
  TEST_REAL_SIMILAR(features[11].getMetaValue("inverse_avgFWHM"), 3.78554915619634)
  TEST_REAL_SIMILAR(features[11].getMetaValue("avgSNR"), 8.18435900228459)
  TEST_REAL_SIMILAR(features[11].getMetaValue("avgFWHM"), 0.264162465929985)

  TEST_EQUAL(features[15].getMetaValue("transition_name"), "glu-L.glu-L_m1-1")
  TEST_REAL_SIMILAR(features[15].getIntensity(), 13.5799)
  TEST_REAL_SIMILAR(features[15].getMetaValue("log10_total_tic"), 5.49089446225569)
  TEST_REAL_SIMILAR(features[15].getMetaValue("inverse_avgFWHM"), 3.53583924309525)
  TEST_REAL_SIMILAR(features[15].getMetaValue("avgSNR"), 4.55314284068408)
  TEST_REAL_SIMILAR(features[15].getMetaValue("avgFWHM"), 0.282818287611008)

  TEST_EQUAL(features[19].getMetaValue("transition_name"), "skm.skm_m4-4")
  TEST_REAL_SIMILAR(features[19].getIntensity(), 10.5746)
  TEST_REAL_SIMILAR(features[19].getMetaValue("log10_total_tic"), 6.60354130105922)
  TEST_REAL_SIMILAR(features[19].getMetaValue("inverse_avgFWHM"), 2.02868912178847)
  TEST_REAL_SIMILAR(features[19].getMetaValue("avgSNR"), 1.94235549504842)
  TEST_REAL_SIMILAR(features[19].getMetaValue("avgFWHM"), 0.492929147822516)

  features.pop_back();
  TEST_EXCEPTION(Exception::InvalidSize, tse.scoreSpectra(annotated_spectra, picked_spectra, features, scored_spectra));
}
END_SECTION

START_SECTION(void scoreSpectra(
  const std::vector<MSSpectrum>& annotated_spectra,
  const std::vector<MSSpectrum>& picked_spectra,
  std::vector<MSSpectrum>& scored_spectra
) const)
{
  TargetedSpectraExtractor tse;
  Param params = tse.getParameters();
  params.setValue("GaussFilter:gaussian_width", 0.25);
  params.setValue("peak_height_min", 15000.0);
  params.setValue("peak_height_max", 110000.0);
  params.setValue("fwhm_threshold", 0.23);
  tse.setParameters(params);

  vector<MSSpectrum> annotated_spectra;
  const vector<MSSpectrum>& spectra = experiment.getSpectra();

  tse.annotateSpectra(spectra, targeted_exp, annotated_spectra);

  vector<MSSpectrum> picked_spectra(annotated_spectra.size());
  for (Size i = 0; i < annotated_spectra.size(); ++i)
  {
    tse.pickSpectrum(annotated_spectra[i], picked_spectra[i]);
  }

  for (Int i = annotated_spectra.size() - 1; i >= 0; --i)
  {
    if (picked_spectra[i].empty())
    {
      annotated_spectra.erase(annotated_spectra.begin() + i);
      picked_spectra.erase(picked_spectra.begin() + i);
    }
  }
  TEST_EQUAL(annotated_spectra.size(), 20)
  TEST_EQUAL(annotated_spectra.size(), picked_spectra.size())

  vector<MSSpectrum> scored_spectra;
  tse.scoreSpectra(annotated_spectra, picked_spectra, scored_spectra);

  TEST_EQUAL(scored_spectra.size(), 20)
  TEST_EQUAL(scored_spectra.size(), annotated_spectra.size())

  TEST_EQUAL(scored_spectra[0].getName(), "met-L.met-L_m0-0")
  TEST_REAL_SIMILAR(scored_spectra[0].getFloatDataArrays()[1][0], 15.2046270370483) // score
  TEST_REAL_SIMILAR(scored_spectra[0].getFloatDataArrays()[2][0], 5.3508939743042)  // total tic
  TEST_REAL_SIMILAR(scored_spectra[0].getFloatDataArrays()[3][0], 3.96267318725586) // inverse average fwhm
  TEST_REAL_SIMILAR(scored_spectra[0].getFloatDataArrays()[4][0], 5.89106035232544) // average snr

  TEST_EQUAL(scored_spectra[4].getName(), "asp-L.asp-L_m1-0")
  TEST_REAL_SIMILAR(scored_spectra[4].getFloatDataArrays()[1][0], 10.8893)
  TEST_REAL_SIMILAR(scored_spectra[4].getFloatDataArrays()[2][0], 6.49946)
  TEST_REAL_SIMILAR(scored_spectra[4].getFloatDataArrays()[3][0], 2.65215)
  TEST_REAL_SIMILAR(scored_spectra[4].getFloatDataArrays()[4][0], 1.73772)

  TEST_EQUAL(scored_spectra[8].getName(), "asp-L.asp-L_m2-1")
  TEST_REAL_SIMILAR(scored_spectra[8].getFloatDataArrays()[1][0], 16.1929)
  TEST_REAL_SIMILAR(scored_spectra[8].getFloatDataArrays()[2][0], 5.52142)
  TEST_REAL_SIMILAR(scored_spectra[8].getFloatDataArrays()[3][0], 3.44492)
  TEST_REAL_SIMILAR(scored_spectra[8].getFloatDataArrays()[4][0], 7.22662)

  TEST_EQUAL(scored_spectra[11].getName(), "asp-L.asp-L_m2-2")
  TEST_REAL_SIMILAR(scored_spectra[11].getFloatDataArrays()[1][0], 17.4552)
  TEST_REAL_SIMILAR(scored_spectra[11].getFloatDataArrays()[2][0], 5.48532)
  TEST_REAL_SIMILAR(scored_spectra[11].getFloatDataArrays()[3][0], 3.78555)
  TEST_REAL_SIMILAR(scored_spectra[11].getFloatDataArrays()[4][0], 8.18436)

  TEST_EQUAL(scored_spectra[15].getName(), "glu-L.glu-L_m1-1")
  TEST_REAL_SIMILAR(scored_spectra[15].getFloatDataArrays()[1][0], 13.5799)
  TEST_REAL_SIMILAR(scored_spectra[15].getFloatDataArrays()[2][0], 5.49089)
  TEST_REAL_SIMILAR(scored_spectra[15].getFloatDataArrays()[3][0], 3.53584)
  TEST_REAL_SIMILAR(scored_spectra[15].getFloatDataArrays()[4][0], 4.55314)

  TEST_EQUAL(scored_spectra[19].getName(), "skm.skm_m4-4")
  TEST_REAL_SIMILAR(scored_spectra[19].getFloatDataArrays()[1][0], 10.5746)
  TEST_REAL_SIMILAR(scored_spectra[19].getFloatDataArrays()[2][0], 6.60354)
  TEST_REAL_SIMILAR(scored_spectra[19].getFloatDataArrays()[3][0], 2.02869)
  TEST_REAL_SIMILAR(scored_spectra[19].getFloatDataArrays()[4][0], 1.94236)
}
END_SECTION

START_SECTION(void selectSpectra(
  const std::vector<MSSpectrum>& scored_spectra,
  const FeatureMap& features,
  std::vector<MSSpectrum>& selected_spectra,
  FeatureMap& selected_features,
  const bool compute_features = true
) const)
{
  const double min_select_score = 15.0;
  TargetedSpectraExtractor tse;
  Param params = tse.getParameters();
  params.setValue("min_select_score", min_select_score);
  params.setValue("GaussFilter:gaussian_width", 0.25);
  params.setValue("peak_height_min", 15000.0);
  params.setValue("peak_height_max", 110000.0);
  params.setValue("fwhm_threshold", 0.23);
  tse.setParameters(params);

  const std::vector<MSSpectrum>& spectra = experiment.getSpectra();
  std::vector<MSSpectrum> annotated;
  FeatureMap features;
  tse.annotateSpectra(spectra, targeted_exp, annotated, features);
  std::vector<MSSpectrum> picked(annotated.size());
  for (Size i = 0; i < annotated.size(); ++i)
  {
    tse.pickSpectrum(annotated[i], picked[i]);
  }
  for (Int i = annotated.size() - 1; i >= 0; --i)
  {
    if (picked[i].empty())
    {
      annotated.erase(annotated.begin() + i);
      picked.erase(picked.begin() + i);
      features.erase(features.begin() + i);
    }
  }
  std::vector<MSSpectrum> scored;
  tse.scoreSpectra(annotated, picked, features, scored);

  std::vector<MSSpectrum> selected_spectra;
  FeatureMap selected_features;

  tse.selectSpectra(scored, features, selected_spectra, selected_features);
  TEST_EQUAL(selected_spectra.size(), 3)
  TEST_EQUAL(selected_spectra.size(), selected_features.size())
  for (Size i = 0; i < selected_spectra.size(); ++i)
  {
    TEST_NOT_EQUAL(selected_spectra[i].getName(), "")
    TEST_EQUAL(selected_spectra[i].getName(), selected_features[i].getMetaValue("transition_name"))
    TEST_EQUAL(selected_spectra[i].getFloatDataArrays()[1][0], selected_features[i].getIntensity())
    TEST_EQUAL(selected_spectra[i].getFloatDataArrays()[1][0] >= min_select_score, true)
  }

  vector<MSSpectrum>::const_iterator it;
  it = findSpectrumByName(selected_spectra, "asp-L.asp-L_m2-1");
  TEST_REAL_SIMILAR(it->getFloatDataArrays()[1][0], 17.4552230834961)
  it = findSpectrumByName(selected_spectra, "met-L.met-L_m0-0");
  TEST_REAL_SIMILAR(it->getFloatDataArrays()[1][0], 16.0294418334961)
  it = findSpectrumByName(selected_spectra, "asp-L.asp-L_m2-2");
  TEST_REAL_SIMILAR(it->getFloatDataArrays()[1][0], 17.4552)

  features.pop_back();
  TEST_EXCEPTION(Exception::InvalidSize, tse.selectSpectra(scored, features, selected_spectra, selected_features));
}
END_SECTION

START_SECTION(void selectSpectra(
  const std::vector<MSSpectrum>& scored_spectra,
  std::vector<MSSpectrum>& selected_spectra
) const)
{
  const double min_select_score = 15.0;
  TargetedSpectraExtractor tse;
  Param params = tse.getParameters();
  params.setValue("min_select_score", min_select_score);
  params.setValue("GaussFilter:gaussian_width", 0.25);
  params.setValue("peak_height_min", 15000.0);
  params.setValue("peak_height_max", 110000.0);
  params.setValue("fwhm_threshold", 0.23);
  tse.setParameters(params);

  const std::vector<MSSpectrum>& spectra = experiment.getSpectra();
  std::vector<MSSpectrum> annotated;
  tse.annotateSpectra(spectra, targeted_exp, annotated);
  std::vector<MSSpectrum> picked(annotated.size());
  for (Size i = 0; i < annotated.size(); ++i)
  {
    tse.pickSpectrum(annotated[i], picked[i]);
  }
  for (Int i = annotated.size() - 1; i >= 0; --i)
  {
    if (picked[i].empty())
    {
      annotated.erase(annotated.begin() + i);
      picked.erase(picked.begin() + i);
    }
  }
  std::vector<MSSpectrum> scored;
  tse.scoreSpectra(annotated, picked, scored);

  std::vector<MSSpectrum> selected_spectra;

  tse.selectSpectra(scored, selected_spectra);
  TEST_EQUAL(selected_spectra.size(), 3)
  for (Size i = 0; i < selected_spectra.size(); ++i)
  {
    TEST_NOT_EQUAL(selected_spectra[i].getName(), "")
    TEST_EQUAL(selected_spectra[i].getFloatDataArrays()[1][0] >= min_select_score, true)
  }

  vector<MSSpectrum>::const_iterator it;
  it = findSpectrumByName(selected_spectra, "asp-L.asp-L_m2-1");
  TEST_REAL_SIMILAR(it->getFloatDataArrays()[1][0], 17.4552230834961)
  it = findSpectrumByName(selected_spectra, "met-L.met-L_m0-0");
  TEST_REAL_SIMILAR(it->getFloatDataArrays()[1][0], 16.0294418334961)
  it = findSpectrumByName(selected_spectra, "asp-L.asp-L_m2-2");
  TEST_REAL_SIMILAR(it->getFloatDataArrays()[1][0], 17.4552)
}
END_SECTION

START_SECTION(void extractSpectra(
  const MSExperiment& experiment,
  const TargetedExperiment& targeted_exp,
  std::vector<MSSpectrum>& extracted_spectra,
  FeatureMap& extracted_features,
  const bool compute_features = true
) const)
{
  TargetedSpectraExtractor tse;
  Param params = tse.getParameters();
  params.setValue("min_select_score", 15.0);
  params.setValue("GaussFilter:gaussian_width", 0.25);
  params.setValue("peak_height_min", 15000.0);
  params.setValue("peak_height_max", 110000.0);
  params.setValue("fwhm_threshold", 0.23);
  tse.setParameters(params);

  vector<MSSpectrum> extracted_spectra;
  FeatureMap extracted_features;
  tse.extractSpectra(experiment, targeted_exp, extracted_spectra, extracted_features);

  TEST_EQUAL(extracted_spectra.size(), extracted_features.size())
  TEST_EQUAL(extracted_spectra.size(), 3)

  vector<MSSpectrum>::const_iterator it;
  it = findSpectrumByName(extracted_spectra, "asp-L.asp-L_m2-1");
  TEST_REAL_SIMILAR(it->getFloatDataArrays()[1][0], 17.4552230834961)
  it = findSpectrumByName(extracted_spectra, "met-L.met-L_m0-0");
  TEST_REAL_SIMILAR(it->getFloatDataArrays()[1][0], 16.0294418334961)
  it = findSpectrumByName(extracted_spectra, "asp-L.asp-L_m2-2");
  TEST_REAL_SIMILAR(it->getFloatDataArrays()[1][0], 17.4552)
}
END_SECTION

START_SECTION(void extractSpectra(
  const MSExperiment& experiment,
  const TargetedExperiment& targeted_exp,
  std::vector<MSSpectrum>& extracted_spectra
) const)
{
  TargetedSpectraExtractor tse;
  Param params = tse.getParameters();
  params.setValue("min_select_score", 15.0);
  params.setValue("GaussFilter:gaussian_width", 0.25);
  params.setValue("peak_height_min", 15000.0);
  params.setValue("peak_height_max", 110000.0);
  params.setValue("fwhm_threshold", 0.23);
  tse.setParameters(params);

  vector<MSSpectrum> extracted_spectra;
  tse.extractSpectra(experiment, targeted_exp, extracted_spectra);

  TEST_EQUAL(extracted_spectra.size(), 3)

  vector<MSSpectrum>::const_iterator it;
  it = findSpectrumByName(extracted_spectra, "asp-L.asp-L_m2-1");
  TEST_REAL_SIMILAR(it->getFloatDataArrays()[1][0], 17.4552230834961)
  it = findSpectrumByName(extracted_spectra, "met-L.met-L_m0-0");
  TEST_REAL_SIMILAR(it->getFloatDataArrays()[1][0], 16.0294418334961)
  it = findSpectrumByName(extracted_spectra, "asp-L.asp-L_m2-2");
  TEST_REAL_SIMILAR(it->getFloatDataArrays()[1][0], 17.4552)
}
END_SECTION

START_SECTION(void matchSpectrum(
  const MSSpectrum& input_spectrum,
  const MSExperiment& library,
  Comparator& cmp,
  std::vector<Match>& matches
))
{
  // MS Library offered by: MoNa - MassBank of North America
  // Title: GC-MS Spectra
  // http://mona.fiehnlab.ucdavis.edu/downloads
  // https://creativecommons.org/licenses/by/4.0/legalcode
  // Changes made: Only a very small subset of spectra is reproduced

  const String msp_path = OPENMS_GET_TEST_DATA_PATH("MoNA-export-GC-MS_Spectra_reduced_TSE_matchSpectrum.msp");
  const String gcms_fullscan_path = OPENMS_GET_TEST_DATA_PATH("TargetedSpectraExtractor_matchSpectrum_GCMS.mzML");
  const String target_list_path = OPENMS_GET_TEST_DATA_PATH("TargetedSpectraExtractor_matchSpectrum_traML.csv");
  MzMLFile mzml;
  MSExperiment gcms_experiment;
  TransitionTSVFile tsv_reader;
  TargetedExperiment targeted_exp;
  mzml.load(gcms_fullscan_path, gcms_experiment);
  Param tsv_params = tsv_reader.getParameters();
  tsv_params.setValue("retentionTimeInterpretation", "seconds");
  tsv_reader.setParameters(tsv_params);
  tsv_reader.convertTSVToTargetedExperiment(target_list_path.c_str(), FileTypes::CSV, targeted_exp);
  TargetedSpectraExtractor tse;
  Param params = tse.getParameters();
  params.setValue("rt_window", 2.0);
  params.setValue("min_select_score", 0.1);
  params.setValue("GaussFilter:gaussian_width", 0.1);
  params.setValue("PeakPickerHiRes:signal_to_noise", 0.01);
  params.setValue("top_matches_to_report", 2);
  params.setValue("min_match_score", 0.51);
  tse.setParameters(params);

  TEST_EQUAL(gcms_experiment.getSpectra().size(), 11)

  vector<MSSpectrum> extracted_spectra;
  FeatureMap extracted_features;
  tse.extractSpectra(gcms_experiment, targeted_exp, extracted_spectra, extracted_features);

  TEST_EQUAL(extracted_spectra.size(), 18)

  MSExperiment library;
  MSPGenericFile mse(msp_path, library);

  TEST_EQUAL(library.getSpectra().size(), 21)

  vector<TargetedSpectraExtractor::Match> matches;

  TargetedSpectraExtractor::BinnedSpectrumComparator cmp;
  std::map<String,DataValue> options = {
    {"bin_size", 1.0},
    {"peak_spread", 0.0},
    {"bin_offset", 0.4}
  };
  cmp.init(library.getSpectra(), options);

  tse.matchSpectrum(extracted_spectra[0], cmp, matches);
  TEST_EQUAL(matches.size() >= 2, true)
  TEST_EQUAL(matches[0].score >= matches[1].score, true)

  tse.matchSpectrum(extracted_spectra[4], cmp, matches);
  TEST_EQUAL(matches.size() >= 2, true)
  TEST_EQUAL(matches[0].score >= matches[1].score, true)

  tse.matchSpectrum(extracted_spectra[8], cmp, matches);
  TEST_EQUAL(matches.size() >= 2, true)
  TEST_EQUAL(matches[0].score >= matches[1].score, true)

  tse.matchSpectrum(extracted_spectra[9], cmp, matches);
  TEST_EQUAL(matches.size() >= 2, true)
  TEST_EQUAL(matches[0].score >= matches[1].score, true)

  tse.matchSpectrum(extracted_spectra[13], cmp, matches);
  TEST_EQUAL(matches.size() >= 2, true)
  TEST_EQUAL(matches[0].score >= matches[1].score, true)

  tse.matchSpectrum(extracted_spectra[17], cmp, matches);
  TEST_EQUAL(matches.size() >= 2, true)
  TEST_EQUAL(matches[0].score >= matches[1].score, true)
}
END_SECTION

START_SECTION(void targetedMatching(
  const std::vector<MSSpectrum>& spectra,
  const MSExperiment& library,
  Comparator& cmp,
  FeatureMap& features
))
{
  // MS Library offered by: MoNa - MassBank of North America
  // Title: GC-MS Spectra
  // http://mona.fiehnlab.ucdavis.edu/downloads
  // https://creativecommons.org/licenses/by/4.0/legalcode
  // Changes made: Only a very small subset of spectra is reproduced

  const String msp_path = OPENMS_GET_TEST_DATA_PATH("MoNA-export-GC-MS_Spectra_reduced_TSE_matchSpectrum.msp");
  const String gcms_fullscan_path = OPENMS_GET_TEST_DATA_PATH("TargetedSpectraExtractor_matchSpectrum_GCMS.mzML");
  const String target_list_path = OPENMS_GET_TEST_DATA_PATH("TargetedSpectraExtractor_matchSpectrum_traML.csv");
  MzMLFile mzml;
  MSExperiment gcms_experiment;
  TransitionTSVFile tsv_reader;
  TargetedExperiment targeted_exp;
  mzml.load(gcms_fullscan_path, gcms_experiment);
  Param tsv_params = tsv_reader.getParameters();
  tsv_params.setValue("retentionTimeInterpretation", "seconds");
  tsv_reader.setParameters(tsv_params);
  tsv_reader.convertTSVToTargetedExperiment(target_list_path.c_str(), FileTypes::CSV, targeted_exp);
  TargetedSpectraExtractor tse;
  Param params = tse.getParameters();
  params.setValue("rt_window", 2.0);
  params.setValue("min_select_score", 0.1);
  params.setValue("GaussFilter:gaussian_width", 0.1);
  params.setValue("PeakPickerHiRes:signal_to_noise", 0.01);
  params.setValue("top_matches_to_report", 2);
  params.setValue("min_match_score", 0.51);
  tse.setParameters(params);

  TEST_EQUAL(gcms_experiment.getSpectra().size(), 11)

  vector<MSSpectrum> extracted_spectra;
  FeatureMap extracted_features;
  tse.extractSpectra(gcms_experiment, targeted_exp, extracted_spectra, extracted_features);

  TEST_EQUAL(extracted_spectra.size(), 18)

  MSExperiment library;
  MSPGenericFile mse(msp_path, library);

  TEST_EQUAL(library.getSpectra().size(), 21)

  TargetedSpectraExtractor::BinnedSpectrumComparator cmp;
  std::map<String,DataValue> options = {
    {"bin_size", 1.0},
    {"peak_spread", 0.0},
    {"bin_offset", 0.4}
  };
  cmp.init(library.getSpectra(), options);

  tse.targetedMatching(extracted_spectra, cmp, extracted_features);
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
