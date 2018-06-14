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
// #include <OpenMS/KERNEL/SpectrumHelper.h>
#include <fstream>
#include <regex>

namespace OpenMS
{
  MSPMetaboFile::MSPMetaboFile(const String& filename, MSExperiment& library)
  {
    load(filename, library);
  }

  void MSPMetaboFile::load(const String& filename, MSExperiment& library)
  {
    // TODO: remove times debug info
    std::clock_t start;
    start = std::clock();
    LOG_INFO << "\nLoading spectra from .msp file. Please wait." << std::endl;
    loaded_spectra_names_.clear();
    std::ifstream ifs(filename, std::ifstream::in);
    if (!ifs.is_open())
    {
      throw Exception::FileNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename);
    }
    const Size BUFSIZE { 65536 };
    char line[BUFSIZE];
    library.clear(true);
    MSSpectrum spectrum;
    bool adding_spectrum { false }; // to avoid calling `.addSpectrum()` on empty/invalid spectra

    std::cmatch m;
    std::regex re_name("^Name: (.+)");
    std::regex re_points_line("^\\d");
    std::regex re_point("(\\d+) (\\d+); ");
    std::regex re_metadatum(" *([^;\r\n]+): ([^;\r\n]+)");

    while (!ifs.eof())
    {
      ifs.getline(line, BUFSIZE);
      // Peaks
      if (std::regex_search(line, m, re_points_line))
      {
        // LOG_DEBUG << "re_points_line\n";
        std::regex_search(line, m, re_point);
        do
        {
          const double position { std::stod(m[1]) };
          const double intensity { std::stod(m[2]) };
          spectrum.push_back( Peak1D(position, intensity) );
          // LOG_DEBUG << position << " " << intensity << "; ";
        } while ( std::regex_search(m[0].second, m, re_point) );
      }
      // Name
      else if (std::regex_search(line, m, re_name))
      {
        addSpectrumToLibrary(spectrum, adding_spectrum, library);
        // LOG_DEBUG << "\n\nName: " << m[1] << "\n";
        spectrum.clear(true);
        spectrum.setName( String(m[1]) );
        adding_spectrum = true;
      }
      // Other metadata
      else if (std::regex_search(line, m, re_metadatum))
      {
        pushParsedInfoToNamedDataArray(spectrum, String(m[1]), String(m[2]));
        while (std::regex_search(m[0].second, m, re_metadatum))
        {
          pushParsedInfoToNamedDataArray(spectrum, String(m[1]), String(m[2]));
        }
      }
    }
    // To make sure a spectrum is added even if no empty line is present before EOF
    addSpectrumToLibrary(spectrum, adding_spectrum, library);
    ifs.close();
    LOG_INFO << "Loading spectra from .msp file completed." << std::endl;
    std::cout << "PARSE TIME: " << ((std::clock() - start) / (double)CLOCKS_PER_SEC) << std::endl;
  }

  void MSPMetaboFile::pushParsedInfoToNamedDataArray(
    MSSpectrum& spectrum,
    const String& name,
    const String& info
  ) const
  {
    // // LOG_DEBUG << name << ": " << info << "\n";
    // MSSpectrum::StringDataArrays& SDAs = spectrum.getStringDataArrays();
    // MSSpectrum::StringDataArrays::iterator it = getDataArrayByName(SDAs, name);
    // if (it != SDAs.end()) // DataArray with given name already exists
    // {
    //   it->push_back(info);
    // }
    // else // DataArray with given name does not exist, yet. Create it.
    // {
    //   MSSpectrum::StringDataArray sda;
    //   sda.push_back(info);
    //   sda.setName(name);
    //   SDAs.push_back(sda);
    // }

    // LOG_DEBUG << name << ": " << info << "\n";
    try // DataArray with given name already exists
    {
      MSSpectrum::StringDataArray& sda = spectrum.getStringDataArrayByName(name);
      sda.push_back(info);
    }
    catch (Exception::ElementNotFound& e) // DataArray with given name does not exist, yet. Create it.
    {
      MSSpectrum::StringDataArray sda;
      sda.push_back(info);
      sda.setName(name);
      spectrum.getStringDataArrays().push_back(sda);
    }
  }

  void MSPMetaboFile::addSpectrumToLibrary(
    MSSpectrum& spectrum,
    bool& adding_spectrum,
    MSExperiment& library
  )
  {
    if (!adding_spectrum) return;

    // Check that required metadata (Name, Num Peaks) is present
    // Num Peaks is checked later in the code (when verifying for the number of points parsed)
    if (spectrum.getName().empty())
    {
      throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        "The current spectrum misses the Name information.");
    }

    // Ensure Comments metadatum is present and, if not, set it to an empty string
    try
    {
      spectrum.getStringDataArrayByName("Comments");
    }
    catch (const Exception::ElementNotFound& e)
    {
      pushParsedInfoToNamedDataArray(spectrum, "Comments", "");
    }

    // Check that the spectrum is not a duplicate (i.e. already present in `library`)
    std::pair<std::set<String>::const_iterator,bool> inserted = loaded_spectra_names_.insert(spectrum.getName());

    if (inserted.second) // `true` if the insertion took place
    {
      // Check that all expected points are parsed
      const String& num_peaks { spectrum.getStringDataArrayByName("Num Peaks").front() };
      if (spectrum.size() != std::stoul(num_peaks) )
      {
        throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
          num_peaks,
          "The number of points parsed does not coincide with `Num Peaks`."
        );
      }
      library.addSpectrum(spectrum);
      if (loaded_spectra_names_.size() % 20000 == 0)
      {
        LOG_INFO << "Loaded " << loaded_spectra_names_.size() << " spectra..." << std::endl;
      }
    }
    else
    {
      LOG_INFO << "DUPLICATE: " << spectrum.getName() << std::endl;
    }

    adding_spectrum = false;
  }
}
