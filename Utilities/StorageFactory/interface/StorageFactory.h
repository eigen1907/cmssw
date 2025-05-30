#ifndef STORAGE_FACTORY_STORAGE_FACTORY_H
#define STORAGE_FACTORY_STORAGE_FACTORY_H

#include "Utilities/StorageFactory/interface/StorageMaker.h"
#include "Utilities/StorageFactory/interface/LocalFileSystem.h"
#include "Utilities/StorageFactory/interface/IOTypes.h"
#include "Utilities/StorageFactory/interface/IOFlags.h"

#include <memory>
#include <string>
#include <tuple>

#include "oneapi/tbb/concurrent_unordered_map.h"

namespace edm::storage {
  class Storage;
  class StorageProxyMaker;
  class StorageFactory {
  public:
    enum CacheHint { CACHE_HINT_APPLICATION, CACHE_HINT_STORAGE, CACHE_HINT_LAZY_DOWNLOAD, CACHE_HINT_AUTO_DETECT };

    enum ReadHint { READ_HINT_UNBUFFERED, READ_HINT_READAHEAD, READ_HINT_AUTO };

    static const StorageFactory *get(void);
    static StorageFactory *getToModify(void);

    // in GB
    static double defaultMinTempFree() { return 4.; }
    static std::string defaultTempDir();

    ~StorageFactory(void);

    // implicit copy constructor
    // implicit assignment operator

    void setCacheHint(CacheHint value);
    CacheHint cacheHint(void) const;

    void setReadHint(ReadHint value);
    ReadHint readHint(void) const;

    bool enableAccounting(bool enabled);
    bool accounting(void) const;

    void setTimeout(unsigned int timeout);
    unsigned int timeout(void) const;

    void setDebugLevel(unsigned int level);
    unsigned int debugLevel(void) const;

    void setTempDir(const std::string &s, double minFreeSpace);
    std::string tempDir(void) const;
    std::string tempPath(void) const;
    double tempMinFree(void) const;

    void setStorageProxyMakers(std::vector<std::unique_ptr<StorageProxyMaker>> makers);

    void stagein(const std::string &url) const;
    std::unique_ptr<Storage> open(const std::string &url, const int mode = IOFlags::OpenRead) const;
    bool check(const std::string &url, IOOffset *size = nullptr) const;

  private:
    typedef oneapi::tbb::concurrent_unordered_map<std::string, std::shared_ptr<StorageMaker>> MakerTable;

    StorageFactory(void);
    StorageMaker *getMaker(const std::string &proto) const;
    StorageMaker *getMaker(const std::string &url, std::string &protocol, std::string &rest) const;

    // Returns
    // - Storage 's' possibly wrapped in LocalCacheFile
    // - bool telling if LocalCacheFile is used
    std::tuple<std::unique_ptr<Storage>, bool> wrapNonLocalFile(std::unique_ptr<Storage> s,
                                                                const std::string &proto,
                                                                const std::string &path,
                                                                const int mode) const;

    mutable MakerTable m_makers;
    CacheHint m_cacheHint;
    ReadHint m_readHint;
    bool m_accounting;
    double m_tempfree;
    std::string m_temppath;
    std::string m_tempdir;
    std::string m_unusableDirWarnings;
    unsigned int m_timeout;
    unsigned int m_debugLevel;
    LocalFileSystem m_lfs;
    std::vector<std::unique_ptr<StorageProxyMaker>> m_storageProxyMakers_;
    static StorageFactory s_instance;
  };
}  // namespace edm::storage
#endif  // STORAGE_FACTORY_STORAGE_FACTORY_H
