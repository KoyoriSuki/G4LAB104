#include <iostream>
#include <filesystem>
#include <TChain.h>
#include <TFile.h>
#include <TTree.h>

void mergeResults() {
    namespace fs = std::filesystem;

    TChain chain("TreeOfData");
    int n = 0;

    // 遍历当前目录下的所有文件
    for (const auto& entry : fs::directory_iterator(fs::current_path())) {
        if (entry.path().extension() == ".root") {
            std::string path = entry.path().string();
            // 确保只添加以xx.root结尾的文件
            if (1) {
                chain.Add(path.c_str());
                n++;
            }
        }
    }

    if (n == 0) {
        std::cout << "No ROOT files found with suffix 'xx.root'." << std::endl;
        return;
    }

    TFile mergedFile("mergedResults.root", "RECREATE");

    TTree* mergedTree = chain.CloneTree();

    mergedTree->Write();
    mergedFile.Close();

    std::cout << "Merged " << n << " ROOT files with 'xx.root' suffix." << std::endl;
}

