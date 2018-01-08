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
#include <OpenMS/FORMAT/AbsoluteQuantitationStandardsFile.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(AbsoluteQuantitationStandardsFile, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

AbsoluteQuantitationStandardsFile* ptr = 0;
AbsoluteQuantitationStandardsFile* null_ptr = 0;

const String csv_path = OPENMS_GET_TEST_DATA_PATH("150516_calibration_concentrations.csv");

START_SECTION(AbsoluteQuantitationStandardsFile())
{
  ptr = new AbsoluteQuantitationStandardsFile();
  TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION(~AbsoluteQuantitationStandardsFile())
{
  delete ptr;
}
END_SECTION

ptr = new AbsoluteQuantitationStandardsFile();

START_SECTION(void load(
  const String& filename,
  std::vector<AbsoluteQuantitationStandards::runConcentration>& run_concentrations
) const)
{
  std::vector<AbsoluteQuantitationStandards::runConcentration> runs;
  ptr->load(csv_path, runs);
  TEST_EQUAL(runs.size(), 9272)
  TEST_EQUAL(runs[0].sample_name, "150516_CM1_Level1")
  TEST_EQUAL(runs[0].component_name, "23dpg.23dpg_1.Light")
  TEST_EQUAL(runs[0].IS_component_name, "23dpg.23dpg_1.Heavy")
  TEST_REAL_SIMILAR(runs[0].actual_concentration, 0)
  TEST_REAL_SIMILAR(runs[0].IS_actual_concentration, 1)
  TEST_EQUAL(runs[0].concentration_units, "uM")
  TEST_REAL_SIMILAR(runs[0].dilution_factor, 1)
  TEST_EQUAL(runs[3090].sample_name, "150516_CM2_Level6")
  TEST_EQUAL(runs[3090].component_name, "f6p.f6p_1.Light")
  TEST_EQUAL(runs[3090].IS_component_name, "f6p.f6p_1.Heavy")
  TEST_REAL_SIMILAR(runs[3090].actual_concentration, 0.8)
  TEST_REAL_SIMILAR(runs[3090].IS_actual_concentration, 1)
  TEST_EQUAL(runs[3090].concentration_units, "uM")
  TEST_REAL_SIMILAR(runs[3090].dilution_factor, 1)
  TEST_EQUAL(runs[6180].sample_name, "150516_CM1_Level5")
  TEST_EQUAL(runs[6180].component_name, "damp.damp_2.Light")
  TEST_EQUAL(runs[6180].IS_component_name, "damp.damp_1.Heavy")
  TEST_REAL_SIMILAR(runs[6180].actual_concentration, 0)
  TEST_REAL_SIMILAR(runs[6180].IS_actual_concentration, 1)
  TEST_EQUAL(runs[6180].concentration_units, "uM")
  TEST_REAL_SIMILAR(runs[6180].dilution_factor, 1)
  TEST_EQUAL(runs[9271].sample_name, "150516_CM3_Level9")
  TEST_EQUAL(runs[9271].component_name, "utp.utp_2.Light")
  TEST_EQUAL(runs[9271].IS_component_name, "utp.utp_1.Heavy")
  TEST_REAL_SIMILAR(runs[9271].actual_concentration, 0)
  TEST_REAL_SIMILAR(runs[9271].IS_actual_concentration, 1)
  TEST_EQUAL(runs[9271].concentration_units, "uM")
  TEST_REAL_SIMILAR(runs[9271].dilution_factor, 1)
}
END_SECTION

delete ptr;

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
