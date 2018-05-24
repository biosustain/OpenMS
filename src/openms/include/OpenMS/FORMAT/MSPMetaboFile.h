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

#pragma once

#include <OpenMS/config.h> // OPENMS_DLLAPI
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/KERNEL/MSExperiment.h>

namespace OpenMS
{
  /**
    @brief Load MSP text file and save it into a `MSExperiment`

    An example of the expected format:
    > Name: spectrum_name
    > 35 310; 36 1230; 37 27; 38 303; 47 5240; 
    > 66 203; 67 68; 68 77; 82 63; 83 240; 
    > 136 350; 
  */
  class OPENMS_DLLAPI MSPMetaboFile
  {
public:
    /// Default constructor
    MSPMetaboFile() = default;

    /// Constructor with filename and output experiment
    MSPMetaboFile(const String& filename, MSExperiment& experiment);

    /// Destructor
    ~MSPMetaboFile() = default;

    /// To test private and protected methods
    friend class MSPMetaboFile_friend;

    /**
      @brief Load the file's data and metadata, and save it into a `MSExperiment`.

      @param[in] filename Path to the MSP input file
      @param[out] experiment The variable into which the extracted information will be saved

      @throw FileNotFound is thrown if the file could not be found
    */
    void load(const String& filename, MSExperiment& experiment);

private:
    /**
      TODO: complete docs
    */
    void pushParsedInfoToNamedDataArray(
      MSSpectrum& spectrum,
      const String& name,
      const String& info
    ) const;

    /**
      TODO: complete docs
    */
    void addSpectrumToExperiment(
      MSSpectrum& spectrum,
      bool& adding_spectrum,
      MSExperiment& experiment
    );

    /**
      TODO: complete docs
    */
    const MSSpectrum::StringDataArray& getStringDataArrayByName(
      const MSSpectrum& spectrum,
      const String& name
    ) const;

    /// To keep track of which spectra have already been loaded and avoid duplicates
    std::vector<String> loaded_spectra_names_;
  };

  class MSPMetaboFile_friend
  {
public:
    MSPMetaboFile_friend() = default;
    ~MSPMetaboFile_friend() = default;

    const MSSpectrum::StringDataArray& getStringDataArrayByName(
      const MSSpectrum& spectrum,
      const String& name
    ) const
    {
      return msp_.getStringDataArrayByName(spectrum, name);
    }

    MSPMetaboFile msp_;
  };
}
