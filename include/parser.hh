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

struct Event {
  long nps, eventID;
  long cellID;
  std::vector<double> pos;
  double energy;
  double weight;
  double time;
};

/**
 * @brief class of a single particle's history.
 * consists of multiple events.
 * 
 */
typedef std::list<Event> ParticleHistory;

/**
 * @brief class of all particles' histories in one nps simulation.
 * 
 */
typedef std::vector<ParticleHistory> NPSHistory;

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

struct EventIndices {
  int event, cell, px, py, pz, erg, wt, tme;
  VariableIDNum idNum;
};

class MCNPPTRAC
{
protected:
  long npsRead;
  // Event record;
  NPSHistory npsHistory;

public:
  MCNPPTRAC();
  virtual ~MCNPPTRAC(){};

  /**
     * If the maximum number of histories has not been reached: reads the next
     * history in PTRAC file.
     *
     * @returns true if successful, false otherwise.
     */
  virtual bool readNextNPS(long maxReadNPS) = 0;

  /**
     * Increments the number of histories read so far.
     */
  void incrementNPSRead();
  long getNPSRead();

  NPSHistory const &getNPSHistory() const;
};

class MCNPPTRACBinary : public MCNPPTRAC
{
protected:
  std::ifstream ptracFile;
  EventIndices indices;

public:
  /**
     * @param[in] ptracPath MCNP ptrac file path.
     */
  MCNPPTRACBinary(std::string const &ptracPath);

  /**
     * If the maximum number of histories has not been reached: reads the history
     * in PTRAC file.
     *
     * @returns true if successful, false otherwise.
     */
  bool readNextNPS(long maxReadNPS);

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

