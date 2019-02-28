#include <string>

class ConnectionInfoPublisher
{
public:
  
  ConnectionInfoPublisher(std::string acceptorName,
                          std::string requesterName,
                          int rank,
                          std::string addressDirectory);

  ConnectionInfoPublisher(std::string acceptorName,
                          std::string requesterName,
                          std::string addressDirectory);

protected:

  int const rank = -1;
  std::string const acceptorName;
  std::string const requesterName;
  std::string const addressDirectory;

  /// Returns the file name for the connection information.
  /**
   * It has the form addressDirectory/.precice/first two letters from hash of 
   * (acceptorName, requesterName, rank)/rest of hash.
   */
  std::string getFilename();  
};

/// Reads the connection info for the given participant/rank information
class ConnectionInfoReader : public ConnectionInfoPublisher
{
public:
  using ConnectionInfoPublisher::ConnectionInfoPublisher; // to inherit the constructor

  /// Reads the info from the connection info file. Will block, if the the file is not present.
  std::string read();
};

/// Writes the connection info for the given participant/rank information.
/**
 * The file is removed, when the object is destroyed.
 */
class ConnectionInfoWriter : public ConnectionInfoPublisher
{
public:
  using ConnectionInfoPublisher::ConnectionInfoPublisher; // to inherit the constructor

  /// Removes the connection info file.
  ~ConnectionInfoWriter();

  /// Write the string info, e.g. IP:port to the connection info file
  /**
   * which is determined by acceptorName, requesterName, rank, addressDirectory
   * set at construction.
   */
  void write(std::string const & info);
};
