inline void print_help_msg()
{    
    std::printf("\n    ");
    std::printf("\n    BST_DUNEURO_MEEG");
    std::printf("\n    This script computes lead fields for EEG, MEG and combined MEG/EEG forward problem.");
    std::printf("\n    Based on duneuro softweare (www.duneuro.org)");
    std::printf("\n    Designed to work with Brainstorm Toolbox (https://neuroimage.usc.edu/brainstorm)");
    std::printf("      ");
    std::printf("\n    Created initially by the Duneuro team.");
    std::printf("\n    Adapted and modified by Takfarinas Medani and Juan Garcia-Prieto");
    std::printf("\n    ");
    std::printf("\n    Usage: bst_duneuro_meeg config_file.ini --mode");
    std::printf("\n    ");
    std::printf("\n        - config_file.ini");
    std::printf("\n          Configuration file containing the parmeters of the forward model to compute.");
    std::printf("\n    ");
    std::printf("\n        - mode: {eeg, meg, meeg}");
    std::printf("\n    ");
    std::printf("\n    This application computes the MEG, MEG and combined MEG/EEG transfer matrix and the final leadfields solution");
	std::printf("\n    If the source/sensor models are not modified, the application will search for a previously computed transfer");
    std::printf("\n    matrix (stored in eeg_transfer.dat or meg_transfer.dat binary files.");
    std::printf("\n    The final output leadfiedl martix will be stored in a binary and a text file named eeg_lf.dat/txt and/or");
    std::printf("\n    meg_lf.dat/txt");
    std::printf("\n    ");
    std::printf("\n    ");
    std::printf("\n    Examples:");
    std::printf("\n    ");
    std::printf("\n        - bst_duneuro --help                                 Will print this help.");
    std::printf("\n    ");
    std::printf("\n        - bst_duneuro eeg_config_file.mini --eeg             Will compute solution for eeg.");
    std::printf("\n    ");
    std::printf("\n        - bst_duneuro meg_config_file.mini --meg             Will compute solution for meg.");
    std::printf("\n    ");
    std::printf("\n        - bst_duneuro meeg_config_file.mini --meeg           Will compute solution for both/combined meg and eeg.");
    std::printf("\n    ");
    std::printf("\n    ");
    std::printf("\n=============================================================================");
    std::printf("\n  This function is part of the Brainstorm software:");
    std::printf("\n  https://neuroimage.usc.edu/brainstorm");
    std::printf("\n  ");
    std::printf("\n  Copyright (c)2000-2019 University of Southern California & McGill University");
    std::printf("\n  This software is distributed under the terms of the GNU General Public License");
    std::printf("\n  as published by the Free Software Foundation. Further details on the GPLv3");
    std::printf("\n  license can be found at http://www.gnu.org/copyleft/gpl.html.");
    std::printf("\n  ");
    std::printf("\n  FOR RESEARCH PURPOSES ONLY. THE SOFTWARE IS PROVIDED \"AS IS,\" AND THE");
    std::printf("\n  UNIVERSITY OF SOUTHERN CALIFORNIA AND ITS COLLABORATORS DO NOT MAKE ANY");
    std::printf("\n  WARRANTY, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO WARRANTIES OF");
    std::printf("\n  MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE, NOR DO THEY ASSUME ANY");
    std::printf("\n  LIABILITY OR RESPONSIBILITY FOR THE USE OF THIS SOFTWARE.");
    std::printf("\n  For more information type \"brainstorm license\" at command prompt.");
    std::printf("\n=============================================================================");
    std::printf("\n    ");
}

inline int powInt(int x, int y)
{
    for (int i = 0; i < y; i++)
    {
        x *= 10;
    }
    return x;
}

inline size_t parseInt(std::vector<int> &v)
{
    size_t sum{0};
    int len{v.size()};
    for (int i; i < len; i++)
    {
        int n{v.at(len - (i + 1)) - '0'};
        sum += powInt(n, i);
    }
    return sum;
}

inline bool fileExists(const std::string &fname)
{
    if (FILE *file = fopen(fname.c_str(), "r"))
    {
        fclose(file);
        return true;
    }
    else
    {
        return false;
    }
}

void saveTransferMatrix(const std::string &fname,
                        const std::shared_ptr<duneuro::DenseMatrix<double>> tMat)
{
    std::ofstream fout(fname, std::ios::binary);
    fout << "::" << tMat->rows() << "::" << tMat->cols() << "::";
    fout.write(reinterpret_cast<char *>(tMat->data()),
               tMat->rows() * tMat->cols() * sizeof(tMat->data()[0]));
    fout.close();
}

std::shared_ptr<duneuro::DenseMatrix<double>> readTransferMatrix(const std::string &fname)
{
    std::ifstream fin{fname, std::ios::binary};

    //read header
    char separatorChar = ':';
    char h1, h2;
    fin.read(&h1, sizeof(h1));
    fin.read(&h2, sizeof(h2));
    if (!(h1 == separatorChar && h2 == separatorChar))
    {
        printf("\nSomething went wrong while reading the binary file.\n");
        return 0;
    }
    std::vector<int> nRowsVec;
    std::vector<int> nColsVec;
    char digit;
    fin.read(&digit, sizeof(digit));
    while (digit != separatorChar)
    {
        nRowsVec.push_back(static_cast<int>(digit));
        fin.read(&digit, sizeof(digit));
    }
    fin.read(&h2, 1);
    if (h2 != separatorChar)
    {
        printf("\nSomething went wrong while reading the binary file.\n");
        return 0;
    }
    fin.read(&digit, sizeof(digit));
    while (digit != separatorChar)
    {
        nColsVec.push_back(static_cast<int>(digit));
        fin.read(&digit, sizeof(digit));
    }
    fin.read(&h1, 1);
    size_t nRows{parseInt(nRowsVec)};
    size_t nCols{parseInt(nColsVec)};
    //printf("\nnumber of rows: %d\n", nRows);
    //printf("\nnumber of cols: %d\n", nCols);

    //now read the data
    std::shared_ptr<duneuro::DenseMatrix<double>>
    matrixTransfer(new duneuro::DenseMatrix<double>(nRows, nCols, 0.));
    fin.read(reinterpret_cast<char *>(matrixTransfer->data()), nRows * nCols * sizeof(matrixTransfer->data()[0]));
    fin.close();
    return matrixTransfer;
}

void saveLFfiles(const std::string &fname,
            std::vector<std::vector<double>> &num_transfer,
            const int numDipoles, const int numElectrodes)
{
    //save leadfields as binary and text files
    std::ofstream fidBin(fname + ".dat", std::ios::binary);
    fidBin << "::" << numDipoles << "::" << numElectrodes << "::";
    for (int i = 0; i < numDipoles; ++i)
    {
        fidBin.write(reinterpret_cast<char *>(&num_transfer[i][0]), numElectrodes * sizeof(num_transfer[i][0]));
    }
    fidBin.close();

    std::ofstream fidTxt{fname + ".txt"};
    for (unsigned int i = 0; i < numDipoles; ++i)
    {
        for (unsigned int j = 0; j < numElectrodes; ++j)
            fidTxt << num_transfer[i][j] << " ";
        fidTxt << "\n";
    }
    fidTxt.close();
}
