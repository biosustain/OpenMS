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

#include <OpenMS/FORMAT/TseMSPFile.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/KERNEL/SpectrumHelper.h>
#include <fstream>
#include <regex>

namespace OpenMS
{
  TseMSPFile::TseMSPFile(const String& filename, MSExperiment& experiment)
  {
    load(filename, experiment);
  }

  void TseMSPFile::load(const String& filename, MSExperiment& experiment) const
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
    std::regex re_num_peaks("^Num Peaks: (.+)");
    std::regex re_points_line("^(?:\\d+ \\d+; ?)+");
    std::regex re_point("(\\d+) (\\d+); ");

    // To keep track of which spectra have already been loaded and avoid duplicates
    std::vector<String> loaded_spectra_names;

    std::function<void(void)> add_spectrum_to_experiment =
    [&experiment, &spectrum, &adding_spectrum, &loaded_spectra_names] ()
    {
      const bool loaded_in_library = std::find(
        loaded_spectra_names.cbegin(),
        loaded_spectra_names.cend(),
        spectrum.getName()
      ) != loaded_spectra_names.cend();
      if (adding_spectrum && !loaded_in_library)
      {
        const String& num_peaks { spectrum.getStringDataArrayByName("Num Peaks").front() };
        if (spectrum.size() != std::stoul(num_peaks) )
        {
          throw Exception::ParseError(
            __FILE__,
            __LINE__,
            OPENMS_PRETTY_FUNCTION,
            num_peaks,
            "Not all peaks could be parsed."
          );
        }
        experiment.addSpectrum(spectrum);
        adding_spectrum = false;
        loaded_spectra_names.push_back(spectrum.getName());
      }
    };

    while (!ifs.eof())
    {
      ifs.getline(line, BUFSIZE);
      // TODO: Should metadata all be saved as strings?
      if (std::regex_search(line, m, re_name))
      {
        LOG_DEBUG << std::endl << std::endl << "Name: " << m[1] << std::endl;
        spectrum.clear(true);
        spectrum.setName( String(m[1]) );
        adding_spectrum = true;
      }
      else if (std::regex_search(line, m, re_synon))
      {
        pushParsedInfoToNamedDataArray(spectrum, "Synon", String(m[1]));
      }
      else if (std::regex_search(line, m, re_formula))
      {
        pushParsedInfoToNamedDataArray(spectrum, "Formula", String(m[1]));
      }
      else if (std::regex_search(line, m, re_mw))
      {
        pushParsedInfoToNamedDataArray(spectrum, "MW", String(m[1]));
      }
      else if (std::regex_search(line, m, re_cas_nist))
      {
        pushParsedInfoToNamedDataArray(spectrum, "CAS#", String(m[1]));
        pushParsedInfoToNamedDataArray(spectrum, "NIST#", String(m[2]));
      }
      else if (std::regex_search(line, m, re_db))
      {
        pushParsedInfoToNamedDataArray(spectrum, "DB#", String(m[1]));
      }
      else if (std::regex_search(line, m, re_comments))
      {
        pushParsedInfoToNamedDataArray(spectrum, "Comments", String(m[1]));
      }
      else if (std::regex_search(line, m, re_num_peaks))
      {
        pushParsedInfoToNamedDataArray(spectrum, "Num Peaks", String(m[1]));
      }
      else if (std::regex_search(line, m, re_points_line))
      {
        LOG_DEBUG << "re_points_line" << std::endl;
        std::regex_search(line, m, re_point);
        do
        {
          const double position { std::stod(m[1]) };
          const double intensity { std::stod(m[2]) };
          spectrum.push_back( Peak1D(position, intensity) );
          LOG_DEBUG << position << " " << intensity << "; ";
        } while ( std::regex_search(m[0].second, m, re_point) );
      }
      else if (line[0] == '\r' || line[0] == '\n')
      {
        LOG_DEBUG << std::endl << "empty_line" << std::endl;
        add_spectrum_to_experiment();
      }
    }
    // To make sure a spectrum is added even if no empty line is present before EOF
    add_spectrum_to_experiment();
    ifs.close();
  }

  void TseMSPFile::pushParsedInfoToNamedDataArray(MSSpectrum& spectrum, const String& name, const String& info) const
  {
    LOG_DEBUG << name << ": " << info << std::endl;
    MSSpectrum::StringDataArrays& SDAs = spectrum.getStringDataArrays();
    MSSpectrum::StringDataArrays::iterator it = getDataArrayByName(SDAs, name);
    if (it != SDAs.end()) // DataArray with given name already exists
    {
      it->push_back(info);
    }
    else // DataArray with given name does not exist, yet. Create it.
    {
      MSSpectrum::StringDataArray sda;
      sda.push_back(info);
      sda.setName(name);
      SDAs.push_back(sda);
    }
  }
}
