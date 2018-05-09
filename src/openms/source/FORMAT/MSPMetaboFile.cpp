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

#include <OpenMS/FORMAT/MSPMetaboFile.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/KERNEL/SpectrumHelper.h>
#include <fstream>
#include <regex>

namespace OpenMS
{
  MSPMetaboFile::MSPMetaboFile(const String& filename, MSExperiment& experiment)
  {
    load(filename, experiment);
  }

  void MSPMetaboFile::load(const String& filename, MSExperiment& experiment)
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
    std::regex re_comments("^Comments: (.+)");
    std::regex re_num_peaks("^Num Peaks: (.+)");
    std::regex re_points_line("^(?:\\d+ \\d+; ?)+");
    std::regex re_point("(\\d+) (\\d+); ");
    std::regex re_metadatum(" *([^;\r\n]+): ([^;\r\n]+)");

    while (!ifs.eof())
    {
      ifs.getline(line, BUFSIZE);
      if (std::regex_search(line, m, re_name))
      {
        LOG_DEBUG << std::endl << std::endl << "Name: " << m[1] << std::endl;
        spectrum.clear(true);
        spectrum.setName( String(m[1]) );
        adding_spectrum = true;
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
      else if (std::regex_search(line, m, re_metadatum))
      {
        pushParsedInfoToNamedDataArray(spectrum, String(m[1]), String(m[2]));
        while (std::regex_search(m[0].second, m, re_metadatum))
        {
          pushParsedInfoToNamedDataArray(spectrum, String(m[1]), String(m[2]));
        }
      }
      else if (line[0] == '\r' || line[0] == '\n')
      {
        LOG_DEBUG << std::endl << "empty_line" << std::endl;
        addSpectrumToExperiment(spectrum, adding_spectrum, experiment);
      }
    }
    // To make sure a spectrum is added even if no empty line is present before EOF
    addSpectrumToExperiment(spectrum, adding_spectrum, experiment);
    ifs.close();
  }

  void MSPMetaboFile::pushParsedInfoToNamedDataArray(
    MSSpectrum& spectrum,
    const String& name,
    const String& info
  ) const
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

  void MSPMetaboFile::addSpectrumToExperiment(
    MSSpectrum& spectrum,
    bool& adding_spectrum,
    MSExperiment& experiment
  )
  {
    if (!adding_spectrum)
    {
      return;
    }

    // Check that required metadata (Name, Num Peaks) is present
    // Num Peaks is checked later in the code (when verifying for the number of points parsed)
    if (spectrum.getName().empty())
    {
      throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        "The current spectrum misses the Name information.");
    }

    // Ensure Comments metadatum is present and, if not, set it to empty string
    try
    {
      getStringDataArrayByName(spectrum, "Comments");
    }
    catch (const Exception::ElementNotFound& e)
    {
      pushParsedInfoToNamedDataArray(spectrum, "Comments", "");
    }

    // Check that the spectrum is not a duplicate (i.e. already present in `experiment`)
    std::vector<String>& l = loaded_spectra_names_;
    const bool is_duplicate = std::find(l.cbegin(), l.cend(), spectrum.getName()) != l.cend();

    if (!is_duplicate)
    {
      // Check that all expected points are parsed
      const String& num_peaks { getStringDataArrayByName(spectrum, "Num Peaks").front() };
      if (spectrum.size() != std::stoul(num_peaks) )
      {
        throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
          num_peaks,
          "The number of points parsed does not coincide with `Num Peaks`."
        );
      }
      experiment.addSpectrum(spectrum);
      loaded_spectra_names_.push_back(spectrum.getName());
    }

    adding_spectrum = false;
  };

  const MSSpectrum::StringDataArray& MSPMetaboFile::getStringDataArrayByName(
    const MSSpectrum& spectrum,
    const String& name
  ) const
  {
    MSSpectrum::StringDataArrays::const_iterator it = std::find_if(
      spectrum.getStringDataArrays().cbegin(),
      spectrum.getStringDataArrays().cend(),
      [&name] (const MSSpectrum::StringDataArray& sda) { return sda.getName() == name; }
    );
    if (it == spectrum.getStringDataArrays().cend())
    {
      throw Exception::ElementNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, name);
    }
    return *it;
  }
}
