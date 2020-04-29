#ifndef READ_DATA_H
#define READ_DATA_H


std::vector<int16_t> readFromOneFile(const std::string name, const int &nEmitter) {
    int elements_num = 2048, time_points_num = 3750;
    int offset = (elements_num / 4) * time_points_num * nEmitter * sizeof(int16_t);
    std::ifstream file(name, std::fstream::binary);
    int lenght = time_points_num * elements_num / 4;
    std::vector<int16_t> result;
    if (file) {
        file.seekg(offset, std::fstream::beg);
        for (int i = 0; i < lenght; i++) {
            int16_t temp = 0;
            file.read((char*)&temp, sizeof(int16_t));
            result.push_back(temp);
        }
        return result;
    }
}

std::vector<std::vector<int>> readEmitter(const std::string directory, const int &nEmitter) {
    std::string file1 = directory + "/decode_data_01.bin";
    std::string file2 = directory + "/decode_data_02.bin";
    std::string file3 = directory + "/decode_data_03.bin";
    std::string file4 = directory + "/decode_data_04.bin";
    std::vector<int16_t> result1, result2, result3, result4;
    result1 = readFromOneFile(file1, nEmitter);
    result2 = readFromOneFile(file2, nEmitter);
    result3 = readFromOneFile(file3, nEmitter);
    result4 = readFromOneFile(file4, nEmitter);
    std::vector<std::vector<int>> result(2048);
    for (int i = 0; i < 2048; i++)
        result[i] = std::vector<int>(3750);
    for (int i = 0; i < 3750; i++) {
        for (int j = 0; j < 512; j ++)
            result[j][i] = ((int)result1[i * 512 + j]);
        for (int j = 512; j < 1024; j ++)
            result[j][i] = ((int)result2[i * 512 + j - 512]);
        for (int j = 1024; j < 1536; j ++)
            result[j][i] = ((int)result3[i * 512 + j - 1024]);
        for (int j = 1536; j < 2048; j ++)
            result[j][i] = ((int)result4[i * 512 + j - 1536]);
    }
    return result;
}

#endif // READ_DATA_H