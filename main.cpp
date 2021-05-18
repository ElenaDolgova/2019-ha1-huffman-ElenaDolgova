#include <string>
#include "huffman.hpp"

enum Arguments {
    verbose = 'v',
    compress = 'c',
    decompress = 'd',
};

int main(int argc, char *argv[]) {
    if (argc == 4) {
        char *input_file(&argv[2][0]);
        char *output_file(&argv[3][0]);
        bool needPrint = false;

        std::string mode(argv[1]);
        switch (mode[1]) {
            case compress:
                compression(input_file, output_file, needPrint);
                break;
            case decompress:
                decompression(input_file, output_file, needPrint);
                break;
            default:
                return 1;
        }
    } else if (argc == 5) {
        char *input_file(&argv[3][0]);
        char *output_file(&argv[4][0]);
        char input_verbose(argv[1][1]);
        bool needPrint = false;

        if (input_verbose == verbose) {
            needPrint = true;
        }

        switch (argv[2][1]) {
            case compress:
                compression(input_file, output_file, needPrint);
                break;
            case decompress:
                decompression(input_file, output_file, needPrint);
                break;
            default:
                return 1;
        }
    }

    return 0;
}
