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

#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/FORMAT/MSPTSEFile.h>
#include <fstream>
#include <regex>

namespace OpenMS
{
  MSPTSEFile::MSPTSEFile(const String& filename, MSExperiment& experiment)
  {
    load(filename, experiment);
  }

  void MSPTSEFile::load(const String& filename, MSExperiment& experiment) const
  {
    std::ifstream ifs(filename, std::ifstream::in);
    if (!ifs.is_open())
    {
      throw Exception::FileNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename);
    }
    const Size BUFSIZE { 65536 };
    char line[BUFSIZE];
    experiment.clear(true);
    MSSpectrum spectrum;
    bool adding_spectrum { false }; // to avoid calling `.addSpectrum()` on empty/invalid spectra
    std::cmatch m;
    std::regex re_name("^Name: (.+)");
    std::regex re_synon("^Synon: (.+)");
    std::regex re_formula("^Formula: (.+)");
    std::regex re_mw("^MW: (.+)");
    std::regex re_cas_nist("^CAS#: (.+);  NIST#: (.+)");
    std::regex re_db("^DB#: (.+)");
    std::regex re_comments("^Comments: (.+)");
    // std::regex re_num_peaks("^Num Peaks: .+");
    std::regex re_points_line("^(?:\\d+ \\d+; ?)+");
    std::regex re_point("(\\d+) (\\d+); ");
    while (!ifs.eof())
    {
      ifs.getline(line, BUFSIZE);
      if (std::regex_search(line, m, re_name))
      {
        LOG_DEBUG << std::endl << std::endl << "re_name ";
        spectrum.clear(true);
        spectrum.setName( std::string(m[1]) );
        LOG_DEBUG << spectrum.getName() << std::endl;
        MSSpectrum::StringDataArrays& SDAs = spectrum.getStringDataArrays();
        SDAs.resize(7);
        // TODO: part of this info could be saved as integer instead of string
        // TODO: this number of metadata entries could vary between spectra libraries,
        //       probably will need to make these headers not hardcoded
        SDAs[0].setName("Synon");
        SDAs[1].setName("Formula");
        SDAs[2].setName("MW");
        SDAs[3].setName("CAS#");
        SDAs[4].setName("NIST#");
        SDAs[5].setName("DB#");
        SDAs[6].setName("Comments");
        adding_spectrum = true;
      }
      else if (std::regex_search(line, m, re_synon))
      {
        LOG_DEBUG << "Synon: " << m[1] << std::endl;
        spectrum.getStringDataArrayByName("Synon").push_back( std::string(m[1]) );
      }
      else if (std::regex_search(line, m, re_formula))
      {
        LOG_DEBUG << "Formula: " << m[1] << std::endl;
        spectrum.getStringDataArrayByName("Formula").push_back( std::string(m[1]) );
      }
      else if (std::regex_search(line, m, re_mw))
      {
        LOG_DEBUG << "MW: " << m[1] << std::endl;
        spectrum.getStringDataArrayByName("MW").push_back( std::string(m[1]) );
      }
      else if (std::regex_search(line, m, re_cas_nist))
      {
        LOG_DEBUG << "CAS#: " << m[1] << ";  NIST#: " << m[2] << std::endl;
        spectrum.getStringDataArrayByName("CAS#").push_back( std::string(m[1]) );
        spectrum.getStringDataArrayByName("NIST#").push_back( std::string(m[2]) );
      }
      else if (std::regex_search(line, m, re_db))
      {
        LOG_DEBUG << "DB#: " << m[1] << std::endl;
        spectrum.getStringDataArrayByName("DB#").push_back( std::string(m[1]) );
      }
      else if (std::regex_search(line, m, re_comments))
      {
        LOG_DEBUG << "Comments: " << m[1] << std::endl;
        spectrum.getStringDataArrayByName("Comments").push_back( std::string(m[1]) );
      }
      else if (std::regex_search(line, m, re_points_line))
      {
        LOG_DEBUG << "re_points_line" << std::endl;
        std::regex_search(line, m, re_point);
        do
        {
          LOG_DEBUG << "re_point ";
          const double position { std::stod(m[1]) };
          const double intensity { std::stod(m[2]) };
          spectrum.push_back( Peak1D(position, intensity) );
          LOG_DEBUG << spectrum.back().getPos() << " " << spectrum.back().getIntensity() << "; ";
        } while ( std::regex_search(m[0].second, m, re_point) );
      }
      else if (line[0] == '\r' || line[0] == '\n')
      {
        LOG_DEBUG << std::endl << "empty_line" << std::endl;
        LOG_DEBUG << "line: " << line << std::endl;
        if (adding_spectrum)
        {
          experiment.addSpectrum(spectrum);
          adding_spectrum = false;
        }
      }
    }
    if (adding_spectrum)
    {
      experiment.addSpectrum(spectrum);
    }
    ifs.close();
  }
}
