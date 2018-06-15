// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2018.
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
#include <OpenMS/FORMAT/MSPMetaboFile.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(MSPMetaboFile, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

MSPMetaboFile* ptr = nullptr;
MSPMetaboFile* null_ptr = nullptr;
const String input_filepath = OPENMS_GET_TEST_DATA_PATH("MSPMetaboFile_input.msp");

START_SECTION(MSPMetaboFile())
{
  ptr = new MSPMetaboFile();
  TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION(~MSPMetaboFile())
{
  delete ptr;
}
END_SECTION

START_SECTION(void load(const String& filename, MSExperiment& experiment) const)
{
  MSPMetaboFile msp;
  MSExperiment experiment;
  msp.load(input_filepath, experiment);
  const vector<MSSpectrum>& spectra = experiment.getSpectra();
  TEST_EQUAL(spectra.size(), 3)

  const MSSpectrum& s1 = spectra[0];
  TEST_EQUAL(s1.size(), 14)
  TEST_EQUAL(s1.getName(), "name1 of first")
  TEST_EQUAL(s1.getStringDataArrayByName("Synon")[0], "name2 of 1st")
  TEST_EQUAL(s1.getStringDataArrayByName("Synon")[1], "name3 of firsttt")
  TEST_EQUAL(s1.getStringDataArrayByName("Formula")[0], "A11B22C333")
  TEST_EQUAL(s1.getStringDataArrayByName("MW")[0], "156")
  TEST_EQUAL(s1.getStringDataArrayByName("CAS#")[0], "0123-45-6")
  TEST_EQUAL(s1.getStringDataArrayByName("NIST#")[0], "654321")
  TEST_EQUAL(s1.getStringDataArrayByName("DB#")[0], "1")
  TEST_EQUAL(s1.getStringDataArrayByName("Comments")[0], "Some comment")
  TEST_EQUAL(s1.getStringDataArrayByName("Num Peaks")[0], "14")
  TEST_EQUAL(s1[0].getPos(), 27)
  TEST_EQUAL(s1[0].getIntensity(), 29)
  TEST_EQUAL(s1[5].getPos(), 60)
  TEST_EQUAL(s1[5].getIntensity(), 41)
  TEST_EQUAL(s1[10].getPos(), 90)
  TEST_EQUAL(s1[10].getIntensity(), 168)
  TEST_EQUAL(s1[13].getPos(), 105)
  TEST_EQUAL(s1[13].getIntensity(), 36)

  const MSSpectrum& s2 = spectra[1];
  TEST_EQUAL(s2.size(), 15)
  TEST_EQUAL(s2.getName(), "name1 of second")
  TEST_EQUAL(s2.getStringDataArrayByName("Synon")[0], "name2 of 2nd")
  TEST_EQUAL(s2.getStringDataArrayByName("Synon")[1], "name3 of seconddd")
  TEST_EQUAL(s2.getStringDataArrayByName("Formula")[0], "A44B55C666")
  TEST_EQUAL(s2.getStringDataArrayByName("MW")[0], "589")
  TEST_EQUAL(s2.getStringDataArrayByName("CAS#")[0], "3210-45-6")
  TEST_EQUAL(s2.getStringDataArrayByName("NIST#")[0], "789564")
  TEST_EQUAL(s2.getStringDataArrayByName("DB#")[0], "2")
  TEST_EQUAL(s2.getStringDataArrayByName("Comments")[0], "Some other comment")
  TEST_EQUAL(s2.getStringDataArrayByName("Num Peaks")[0], "15")
  TEST_EQUAL(s2[0].getPos(), 27)
  TEST_EQUAL(s2[0].getIntensity(), 29)
  TEST_EQUAL(s2[5].getPos(), 260)
  TEST_EQUAL(s2[5].getIntensity(), 41)
  TEST_EQUAL(s2[10].getPos(), 290)
  TEST_EQUAL(s2[10].getIntensity(), 168)
  TEST_EQUAL(s2[14].getPos(), 310)
  TEST_EQUAL(s2[14].getIntensity(), 20)

  const MSSpectrum& s3 = spectra[2];
  TEST_EQUAL(s3.size(), 16)
  TEST_EQUAL(s3.getName(), "name1 of third")
  TEST_EQUAL(s3.getStringDataArrayByName("Synon")[0], "name2 of 3rd")
  TEST_EQUAL(s3.getStringDataArrayByName("Synon")[1], "name3 of thirddd")
  TEST_EQUAL(s3.getStringDataArrayByName("Formula")[0], "A12B12C123")
  TEST_EQUAL(s3.getStringDataArrayByName("MW")[0], "562")
  TEST_EQUAL(s3.getStringDataArrayByName("CAS#")[0], "4210-47-4")
  TEST_EQUAL(s3.getStringDataArrayByName("NIST#")[0], "749514")
  TEST_EQUAL(s3.getStringDataArrayByName("DB#")[0], "3")
  TEST_EQUAL(s3.getStringDataArrayByName("Num Peaks")[0], "16")
  TEST_EQUAL(s3[0].getPos(), 27)
  TEST_EQUAL(s3[0].getIntensity(), 29)
  TEST_EQUAL(s3[5].getPos(), 260)
  TEST_EQUAL(s3[5].getIntensity(), 41)
  TEST_EQUAL(s3[10].getPos(), 290)
  TEST_EQUAL(s3[10].getIntensity(), 168)
  TEST_EQUAL(s3[14].getPos(), 310)
  TEST_EQUAL(s3[14].getIntensity(), 20)
  TEST_EQUAL(s3[15].getPos(), 111)
  TEST_EQUAL(s3[15].getIntensity(), 44)
}
END_SECTION

START_SECTION(void pushParsedInfoToNamedDataArray(
  MSSpectrum& spectrum,
  const String& name,
  const String& info
) const)
{
  MSPMetaboFile_friend msp_f;
  MSSpectrum spectrum;

  const String field_synon { "Synon" };
  const String synon1 { "foo" };
  const String synon2 { "bar" };
  msp_f.pushParsedInfoToNamedDataArray(spectrum, field_synon, synon1);
  TEST_EQUAL(spectrum.getStringDataArrays().size(), 1)
  MSSpectrum::StringDataArray& sda_synon = spectrum.getStringDataArrayByName(field_synon);
  TEST_EQUAL(sda_synon.size(), 1)
  TEST_STRING_EQUAL(sda_synon[0], synon1)
  msp_f.pushParsedInfoToNamedDataArray(spectrum, field_synon, synon2);
  TEST_EQUAL(spectrum.getStringDataArrays().size(), 1)
  TEST_EQUAL(sda_synon.size(), 2)
  TEST_STRING_EQUAL(sda_synon[1], synon2)

  const String field_comments { "Comments" };
  const String comment { "seems to work fine" };
  msp_f.pushParsedInfoToNamedDataArray(spectrum, field_comments, comment);
  TEST_EQUAL(spectrum.getStringDataArrays().size(), 2)
  MSSpectrum::StringDataArray& sda_comments = spectrum.getStringDataArrayByName(field_comments);
  TEST_EQUAL(sda_comments.size(), 1)
  TEST_STRING_EQUAL(sda_comments[0], comment)
}
END_SECTION

START_SECTION(void addSpectrumToLibrary(
  MSSpectrum& spectrum,
  MSExperiment& library
))
{
  MSPMetaboFile_friend msp_f;
  MSExperiment lib;

  MSSpectrum spec;
  spec.setName(""); // empty name
  spec.setMetaValue("is_valid", 1);

  TEST_EXCEPTION(Exception::MissingInformation, msp_f.addSpectrumToLibrary(spec, lib))
  TEST_EQUAL(lib.size(), 0)

  spec.setName("foo"); // Num Peaks still absent!
  TEST_EXCEPTION(Exception::ElementNotFound, msp_f.addSpectrumToLibrary(spec, lib))
  TEST_EQUAL(lib.size(), 0)

  msp_f.pushParsedInfoToNamedDataArray(spec, "Num Peaks", "2");
  // Num Peaks is set but raw data poins have not been added
  TEST_EXCEPTION(Exception::ParseError, msp_f.addSpectrumToLibrary(spec, lib))
  TEST_EQUAL(lib.size(), 0)

  spec.push_back(Peak1D(1.0, 2.0));
  spec.push_back(Peak1D(3.0, 4.0)); // now the spectrum is valid
  msp_f.addSpectrumToLibrary(spec, lib);
  TEST_EQUAL(lib.size(), 1)

  spec.setName("bar");
  spec.setMetaValue("is_valid", 1);
  msp_f.addSpectrumToLibrary(spec, lib);
  TEST_EQUAL(lib.size(), 2)

  spec.setMetaValue("is_valid", 1);
  msp_f.addSpectrumToLibrary(spec, lib); // duplicate, won't be added
  TEST_EQUAL(lib.size(), 2)

  spec.setMetaValue("is_valid", 0);
  spec.setName("not a duplicate");
  msp_f.addSpectrumToLibrary(spec, lib);
  TEST_EQUAL(lib.size(), 2)
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
