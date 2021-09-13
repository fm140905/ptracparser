/**
 * @file Parser.hh
 *
 *
 * @brief Binary PTRAC file parser
 *
 * @author J. Faustin, Ming Fang
 * @version 1.1
 */

#pragma once

#include "PTRACFormat.hh"
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <tuple>
#include <utility>
#include <vector>
#include <list>

struct PTRACRecord {
  long pointID, eventID;
  long cellID;
  std::vector<double> point;
  double energy;
  double weight;
  double time;
};

typedef std::list<PTRACRecord> History;

struct VariableIDNum {
  long nbDataNPS;
  long nbDataSrcLong;
  long nbDataSrcDouble;
  long nbDataBnkLong;
  long nbDataBnkDouble;
  long nbDataSufLong;
  long nbDataSufDouble;
  long nbDataColLong;
  long nbDataColDouble;
  long nbDataTerLong;
  long nbDataTerDouble;
};

struct PTRACRecordIndices {
  int event, cell, px, py, pz, erg, wt, tme;
  VariableIDNum idNum;
};

class MCNPPTRAC
{
protected:
  long nbPointsRead;
  // PTRACRecord record;
  History history;

public:
  MCNPPTRAC();
  virtual ~MCNPPTRAC(){};

  /**
     * If the maximum number of read points has not been reached: reads the next
     * particle, event, volume, material, position in PTRAC file.
     *
     * @returns true if successful, false otherwise.
     */
  virtual bool readNextPtracData(long maxReadPoint) = 0;

  /**
     * Increments the number of points read so far.
     */
  void incrementNbPointsRead();
  long getNbPointsRead();

  History const &getPTRACRecord() const;
};

// class MCNPPTRACASCII : public MCNPPTRAC
// {
// protected:
//   std::string currentLine;
//   int nbDataCellMaterialLine;
//   std::ifstream ptracFile;

// public:
//   /**
//      * @param[in] ptracPath MCNP ptrac file path.
//      */
//   MCNPPTRACASCII(std::string const &ptracPath);

//   /**
//      * If the maximum number of read points has not been reached: reads the next
//      * particle, event, volume, material, position in PTRAC file.
//      *
//      * @returns true if successful, false otherwise.
//      */
//   bool readNextPtracData(long maxReadPoint);

// protected:
//   /**
//      * Reads the header lines. Sets the current line at the last header line of PTRAC file.
//      *
//      * @param[in] nHeaderLines The number of header lines in the PTRAC file.
//      */
//   void goThroughHeaderPTRAC(int nHeaderLines);

//   /**
//      * Gets the number of data expected on 2nd line of each particle event data block.
//      * If the PTRAC has the correct format, this information is found in the 5th
//      * header line.
//      *
//      * @param[in] line5 The string containing the data of the 5th header line.
//      * @return the number of integer data stored on the 2nd of each particle event
//      * data block.
//      */
//   int getDataFromLine5Ptrac(const std::string &line5);

//   /**
//      * Checks that the 2nd line of each particle event data block does contain the
//      * cell ID and the material ID. These information are respectively identified as
//      * 17 and 18 by the PTRAC writer. The program exits if it is not the case.
//      *
//      * @param[in] line6 The string containing the data of the 6th header line.
//      * @param[in] nbData The number of data expected on 2nd line of each particle
//      * event data block (given by getDataFromLine5Ptrac function).
//      */
//   void checkDataFromLine6Ptrac(const std::string &line6, int nbData);

//   /**
//      * Reads the point ID number and the event ID number.
//      *
//      * @return A pair containing the point ID and event ID.
//      */
//   std::pair<int, int> readPointEvent();

//   /**
//      * Reads the cell ID number and the material ID number.
//      *
//      * @return a pair containing the volume ID and material ID.
//      */
//   std::pair<int, int> readCellMaterial();

//   /**
//      * Reads the point coordinates (x,y,z).
//      *
//      * @return the point coordinates as a vector of 3 doubles.
//      */
//   std::vector<double> readPoint();

//   /**
//      * Stores the next line in the PTRAC file in the currentLine variable.
//      *
//      *
//      */
//   void getNextLinePtrac();
// };

class MCNPPTRACBinary : public MCNPPTRAC
{
protected:
  std::ifstream ptracFile;
  PTRACRecordIndices indices;

public:
  /**
     * @param[in] ptracPath MCNP ptrac file path.
     */
  MCNPPTRACBinary(std::string const &ptracPath);

  /**
     * If the maximum number of read points has not been reached: reads the next
     * particle, event, volume, material, position in PTRAC file.
     *
     * @returns true if successful, false otherwise.
     */
  bool readNextPtracData(long maxReadPoint);

protected:
  /**
   * Reads the header
   */
  void parseHeader();

  void skipHeader();

  void skipPtracInputData();

  /// return the indices of the field IDs
  void parseVariableIDs();

  void parsePTRACRecord();
};

/*******************************************************
*  utility functions for parsing FORTRAN binary files  *
********************************************************/

template <typename T>
T readBinary(std::istream &stream)
{
  T value;
  stream.read(reinterpret_cast<char *>(&value), sizeof(T) / sizeof(char));
  if (!stream) {
    throw std::logic_error("bad stream state after read");
  }
  return value;
}

inline std::string readBinary(std::istream &file, size_t length)
{
  std::string buffer(length, '\0');
  file.read(&buffer[0], length);
  if (!file) {
    throw std::logic_error("bad stream state after read");
  }
  return buffer;
}

inline std::string readRecord(std::istream &file)
{
  const int rec_len_start = readBinary<int>(file);
  std::string buffer = readBinary(file, rec_len_start);
  const int rec_len_end = readBinary<int>(file);
  if (rec_len_start != rec_len_end) {
    throw std::logic_error("mismatched record length");
  }
  return buffer;
}

inline std::tuple<> reinterpretBuffer(std::istream &, std::tuple<> const &)
{
  return std::make_tuple();
}

template <typename T, typename... Ts>
std::tuple<T, Ts...> reinterpretBuffer(std::istream &stream, std::tuple<T, Ts...> const &)
{
  T value = readBinary<T>(stream);
  return std::tuple_cat(std::make_tuple(value), reinterpretBuffer(stream, std::tuple<Ts...>()));
}

template <typename... Ts>
std::tuple<Ts...> reinterpretBuffer(std::istream &stream)
{
  return reinterpretBuffer<Ts...>(stream, std::tuple<Ts...>());
}

template <typename... Ts>
std::tuple<Ts...> reinterpretBuffer(std::string const &buffer)
{
  std::stringstream stream(buffer);
  return reinterpretBuffer<Ts...>(stream);
}

