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
  // Verifying that Comments is present even if absent within the file
  TEST_EQUAL(s3.getStringDataArrayByName("Comments")[0], "")
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

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
