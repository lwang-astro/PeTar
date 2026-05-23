#pragma once

#include <cstdio>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

class OutputCommitManager {
public:
    enum class CommitMode {
        Append,
        Rename
    };

    struct TmpFileRecord {
        std::string kind;
        int rank;
        std::string tmp_path;
        std::string final_path;
        CommitMode mode;
    };

private:
    bool window_active_;
    std::string window_id_;
    std::vector<TmpFileRecord> records_;

public:
    OutputCommitManager(): window_active_(false), window_id_(), records_() {}

    void beginWindow(const std::string& window_id) {
        window_active_ = true;
        window_id_ = window_id;
        records_.clear();
    }

    void registerTmp(const std::string& kind,
                     const int rank,
                     const std::string& tmp_path,
                     const std::string& final_path,
                     const CommitMode mode) {
        records_.push_back({kind, rank, tmp_path, final_path, mode});
    }

    bool commitWindow() {
        for (const auto& record : records_) {
            bool ok = false;
            switch (record.mode) {
            case CommitMode::Append:
                ok = appendFile(record.tmp_path, record.final_path);
                if (!ok) {
                    std::cerr<<"OutputCommitManager append stage failed: kind="<<record.kind
                             <<" rank="<<record.rank
                             <<" tmp="<<record.tmp_path
                             <<" final="<<record.final_path
                             <<std::endl;
                }
                if (ok) {
                    const bool rm_ok = removeFileSafe(record.tmp_path);
                    if (!rm_ok) {
                        std::cerr<<"OutputCommitManager remove tmp failed: kind="<<record.kind
                                 <<" rank="<<record.rank
                                 <<" tmp="<<record.tmp_path
                                 <<std::endl;
                    }
                    ok = rm_ok;
                }
                break;
            case CommitMode::Rename:
                ok = atomicRename(record.tmp_path, record.final_path);
                if (!ok) {
                    std::cerr<<"OutputCommitManager rename stage failed: kind="<<record.kind
                             <<" rank="<<record.rank
                             <<" tmp="<<record.tmp_path
                             <<" final="<<record.final_path
                             <<std::endl;
                }
                break;
            }
            if (!ok) {
                std::cerr<<"OutputCommitManager commit failure: kind="<<record.kind
                         <<" rank="<<record.rank
                         <<" mode="<<(record.mode==CommitMode::Append?"append":"rename")
                         <<" tmp="<<record.tmp_path
                         <<" final="<<record.final_path
                         <<" tmp_exists="<<(fileExists(record.tmp_path)?1:0)
                         <<" final_exists="<<(fileExists(record.final_path)?1:0)
                         <<std::endl;
                return false;
            }
        }
        records_.clear();
        window_id_.clear();
        window_active_ = false;
        return true;
    }

    void abandonWindow() {
        records_.clear();
        window_id_.clear();
        window_active_ = false;
    }

    bool isWindowActive() const {
        return window_active_;
    }

    const std::string& getWindowId() const {
        return window_id_;
    }

    const std::vector<TmpFileRecord>& getRecords() const {
        return records_;
    }

    static bool appendFile(const std::string& src_tmp, const std::string& dst_final) {
        if (!fileExists(src_tmp)) return true;

        std::FILE* fin = std::fopen(src_tmp.c_str(), "rb");
        if (fin==NULL) {
            std::cerr<<"OutputCommitManager appendFile cannot open src tmp: "<<src_tmp<<std::endl;
            return false;
        }

        std::FILE* fout = std::fopen(dst_final.c_str(), "ab");
        if (fout==NULL) {
            std::cerr<<"OutputCommitManager appendFile cannot open dst final: "<<dst_final<<std::endl;
            std::fclose(fin);
            return false;
        }

        char buffer[1024*1024];
        bool ok = true;
        while (true) {
            const std::size_t n_read = std::fread(buffer, 1, sizeof(buffer), fin);
            if (n_read>0) {
                const std::size_t n_write = std::fwrite(buffer, 1, n_read, fout);
                if (n_write!=n_read) {
                    ok = false;
                    break;
                }
            }
            if (n_read < sizeof(buffer)) {
                if (std::ferror(fin)) ok = false;
                break;
            }
        }

        if (std::fflush(fout)!=0) ok = false;
        std::fclose(fout);
        std::fclose(fin);

        if (!ok) {
            std::cerr<<"OutputCommitManager appendFile write failed: src="<<src_tmp
                     <<" dst="<<dst_final<<std::endl;
        }
        return ok;
    }

    static bool atomicRename(const std::string& src_tmp, const std::string& dst_final) {
        return std::rename(src_tmp.c_str(), dst_final.c_str()) == 0;
    }

    static bool removeFileSafe(const std::string& path) {
        if (!fileExists(path)) return true;
        return std::remove(path.c_str()) == 0;
    }

    static bool fileExists(const std::string& path) {
        std::ifstream fin(path.c_str(), std::ios::binary);
        return fin.good();
    }
};
