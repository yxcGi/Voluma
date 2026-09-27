#pragma once
// 命令行参数：--KEY value 形式，算例用 args.get("KEY", 默认值) 读取（类似 OpenLB 的 fromCLI）。

#include <cstdlib>
#include <iostream>
#include <map>
#include <set>
#include <string>

namespace cfd {

class Args {
public:
    Args(int argc, char** argv) {
        for (int i = 1; i < argc; ++i) {
            std::string k = argv[i];
            if (k.rfind("--", 0) != 0) continue;
            k = k.substr(2);
            if (i + 1 < argc && std::string(argv[i + 1]).rfind("--", 0) != 0)
                kv_[k] = argv[++i];
            else
                kv_[k] = "1";
        }
    }
    bool has(const std::string& k) const { return kv_.count(k) > 0; }
    std::string get(const std::string& k, const std::string& def) const { return use(k) ? kv_.at(k) : def; }
    std::string get(const std::string& k, const char* def) const { return get(k, std::string(def)); }
    double get(const std::string& k, double def) const { return use(k) ? std::atof(kv_.at(k).c_str()) : def; }
    int get(const std::string& k, int def) const { return use(k) ? std::atoi(kv_.at(k).c_str()) : def; }
    bool flag(const std::string& k) const { return use(k) && kv_.at(k) != "0"; }
    // 报告未被读取的参数（多半是拼写错误）
    void warnUnused() const {
        for (auto& [k, v] : kv_)
            if (!used_.count(k)) std::cout << "Warning: unused argument --" << k << '\n';
    }

private:
    bool use(const std::string& k) const {
        used_.insert(k);
        return kv_.count(k) > 0;
    }
    std::map<std::string, std::string> kv_;
    mutable std::set<std::string> used_;
};

} // namespace cfd
