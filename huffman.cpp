#include <set>
#include <map>
#include <vector>
#include <bitset>
#include<climits>
#include <iostream>
#include <fstream>
#include <cmath>
#include <queue>

const u_int32_t MAX_COUNT_BITS_OF_CHAR = 32;
const u_int32_t BYTE_IN_COUNT_OF_CODING_STRING = 4;
const u_int32_t BITS_IN_COUNT_OF_CODING_STRING = 32;
const u_int32_t BYTES_IN_ONE_LINE_OF_TABLE = 1 + 1 + MAX_COUNT_BITS_OF_CHAR / 8;
const u_int32_t BITS_IN_ONE_LINE_OF_TABLE = BYTES_IN_ONE_LINE_OF_TABLE * 8;
const u_long BITS_COUNT = 8;

using namespace std;

namespace {
    class NodeHuffman {
        u_int16_t symbol = 256;
        u_long frequency{};
        NodeHuffman *leftChild = nullptr;
        NodeHuffman *rightChild = nullptr;

    public:
        u_int16_t getSymbol() const {
            return symbol;
        }

        u_long getFrequency() const {
            return frequency;
        }

        NodeHuffman *getLeftChild() const {
            return leftChild;
        }

        NodeHuffman *getRightChild() const {
            return rightChild;
        }

        void setLeftChild(NodeHuffman *child) {
            NodeHuffman::leftChild = child;
        }

        void setRightChild(NodeHuffman *child) {
            NodeHuffman::rightChild = child;
        }

        NodeHuffman() = default;

        NodeHuffman(u_int16_t symbol, long frequency) {
            this->symbol = symbol;
            this->frequency = frequency;
        }

        explicit NodeHuffman(long frequency) {
            this->frequency = frequency;
        }

        ~NodeHuffman() {
            delete leftChild;
            delete rightChild;
        }
    };

    class BTHuffman {
        NodeHuffman *root = nullptr;

    public:
        NodeHuffman *getRoot() const {
            return root;
        }

        BTHuffman() = default;

        BTHuffman(u_int16_t symbol, u_long frequency) : root(new NodeHuffman(symbol, frequency)) {}

        explicit BTHuffman(u_long frequency) : root(new NodeHuffman(frequency)) {}

        void insert(NodeHuffman *nodeHuffman) {
            if (root->getLeftChild() == nullptr) {
                root->setLeftChild(nodeHuffman);
            } else {
                root->setRightChild(nodeHuffman);
            }
        }

        u_long getFrequency() const {
            return root->getFrequency();
        }

        static void insert_node_by_path(NodeHuffman *node, NodeHuffman *newNode, vector<char> *path) {
            for (size_t i = 0; i < path->size() - 1; ++i) {
                if (path->at(i) == '1') {
                    if (node->getRightChild() == nullptr) {
                        node->setRightChild(new NodeHuffman(256, -1));
                        node = node->getRightChild();
                    } else {
                        node = node->getRightChild();
                    }
                } else if (path->at(i) == '0') {
                    if (node->getLeftChild() == nullptr) {
                        node->setLeftChild(new NodeHuffman(256, -1));
                        node = node->getLeftChild();
                    } else {
                        node = node->getLeftChild();
                    }
                }
            }

            if (path->at(path->size() - 1) == '1') {
                node->setRightChild(newNode);
            } else {
                node->setLeftChild(newNode);
            }
        }

        void insert_to_Huffman_by_path(u_int8_t symbol, vector<char> *path) {
            auto *nodeHuffman = new NodeHuffman(symbol, 0);
            if (root == nullptr) {
                root = new NodeHuffman(256, -1);
            }
            insert_node_by_path(root, nodeHuffman, path);
        }

        u_int32_t print_symbol(u_int32_t i, const vector<u_int8_t> *coding_string, vector<u_int8_t> *string) {
            NodeHuffman *node = root;
            while (node != nullptr && node->getSymbol() == 256) {
                if ((*coding_string).at(i) == 1) {
                    i++;
                    node = node->getRightChild();
                    if (node != nullptr && node->getSymbol() != 256) {
                        string->push_back(node->getSymbol());
                        return i;
                    }
                } else {
                    i++;
                    node = node->getLeftChild();
                    if (node != nullptr && node->getSymbol() != 256) {
                        string->push_back(node->getSymbol());
                        return i;
                    }
                }
            }

            return -1;
        }
    };

    void create_alphabet_frequency(const vector<u_int8_t> *buff, map<u_int16_t, u_long> &alphabetWithFrequency) {
        u_long sum = 0;
        for (auto elem: *buff) {
            ++alphabetWithFrequency[elem];
            sum += 1;
        }
        cout << sum << "\n";
    }

    void pre_order(const NodeHuffman *nodeHuffman,
                   map<u_int8_t, pair<bitset<MAX_COUNT_BITS_OF_CHAR>, u_int32_t>> &table,
                   bitset<MAX_COUNT_BITS_OF_CHAR> vv,
                   const u_int8_t bit,
                   const u_int32_t it) {
        if (bit != 2) {
            vv[it] = bit;
        }
        if (nodeHuffman != nullptr) {
            if (nodeHuffman->getSymbol() != 256) {// при инсерте в таблице можно пихнуть u_int16
                table.insert({nodeHuffman->getSymbol(), {vv, it + 1}});
            }
            pre_order(nodeHuffman->getLeftChild(), table, vv, 0, it + 1);
            pre_order(nodeHuffman->getRightChild(), table, vv, 1, it + 1);
        }
    }

    void create_table(const map<u_int16_t, u_long> &alphabetWithFrequency,
                      map<u_int8_t, pair<bitset<MAX_COUNT_BITS_OF_CHAR>, u_int32_t>> &table) {

        auto cmp = [](BTHuffman elem1, BTHuffman elem2) { return elem1.getFrequency() > elem2.getFrequency(); };
        std::priority_queue<BTHuffman, std::vector<BTHuffman>, decltype(cmp)> priorityQueueByFrequency(cmp);

        for (auto alphabet: alphabetWithFrequency) {
            auto btHuffman = new BTHuffman(alphabet.first, alphabet.second);
            priorityQueueByFrequency.push(*btHuffman);
            delete btHuffman;
        }

        while (priorityQueueByFrequency.size() != 1) {
            auto begin = priorityQueueByFrequency.top();
            priorityQueueByFrequency.pop();
            auto nextBegin = priorityQueueByFrequency.top();
            priorityQueueByFrequency.pop();

            auto btHuffmanResult = new BTHuffman(begin.getFrequency() + nextBegin.getFrequency());
            btHuffmanResult->insert(begin.getRoot());
            btHuffmanResult->insert(nextBegin.getRoot());
            priorityQueueByFrequency.push(*btHuffmanResult);
            delete btHuffmanResult;

        }
        bitset<MAX_COUNT_BITS_OF_CHAR> vv;
        pre_order(priorityQueueByFrequency.top().getRoot(), table, vv, 2, -1);
        delete priorityQueueByFrequency.top().getRoot();
    }

    void read_code_char(
            const vector<u_int8_t> &all_bits,
            const u_int8_t count_of_code_char,
            vector<char> *coding_char,
            bitset<MAX_COUNT_BITS_OF_CHAR> &bit,
            const u_int32_t shift,
            const u_int32_t BITS_FIRST_STRING) {
        size_t start = BITS_FIRST_STRING + shift * BITS_IN_ONE_LINE_OF_TABLE + 16;

        for (size_t i = start; i < start + count_of_code_char; ++i) {
            coding_char->push_back('0' + all_bits[i]);
            bit[i - start] = all_bits[i];
        }
    }

    u_int32_t create_compressing_string(
            const vector<u_int8_t> *buff,
            map<u_int8_t, pair<bitset<MAX_COUNT_BITS_OF_CHAR>, u_int32_t>> &table,
            vector<bitset<CHAR_BIT>> &coding_string) {
        u_int32_t countOfWroteBits = 0;
        for (u_int8_t symbol:*buff) {
            bitset<MAX_COUNT_BITS_OF_CHAR> symbol_code = table[symbol].first;
            size_t count_bits_char = table[symbol].second;
            for (size_t i = 0; i < count_bits_char; i++) {
                if (coding_string.size() <= countOfWroteBits / 8) {
                    bitset<CHAR_BIT> newByte;
                    coding_string.push_back(newByte);
                }
                coding_string[countOfWroteBits / 8][7 - countOfWroteBits % 8] = symbol_code[i];
                ++countOfWroteBits;
            }
        }
        cout << coding_string.size() << "\n";

        return countOfWroteBits;
    }


    void write_to_bin(const char *output_file,
                      u_int32_t count_bits_in_string,
                      const vector<bitset<CHAR_BIT>> &code_of_string,
                      const map<u_int8_t, pair<bitset<MAX_COUNT_BITS_OF_CHAR>, u_int32_t>> &table) {
        ofstream out;
        out.open(output_file, ios::binary | ios::out);
        out.write(reinterpret_cast<char *> (&count_bits_in_string),
                  sizeof(count_bits_in_string)); // кол-во бит в стринге

        for (auto &j : code_of_string) {
            unsigned long byte = j.to_ulong();
            out.write(reinterpret_cast<char *> (&byte), 1); // код стринги
        }

        for (auto symbol_code: table) {
            out.write(reinterpret_cast<const char *>(&symbol_code.first), 1); //букву 8 бит

            vector<bitset<CHAR_BIT>> bytesForCodeOfChar;
            u_int8_t count_bits_symbol = symbol_code.second.second;
            out.write(reinterpret_cast<char *> (&count_bits_symbol), 1); //количество бит в кодировке

            bitset<CHAR_BIT> byte;
            for (size_t j = 0; j < MAX_COUNT_BITS_OF_CHAR; j++) {
                byte[7 - j % 8] = symbol_code.second.first[j];
                if ((j + 1) % 8 == 0) {
                    unsigned long code = byte.to_ulong();
                    out.write(reinterpret_cast<char *> (&code), 1); //код Хаффмана буквы
                }
            }
        }
        out.close();
    }

    void create_dictionary(const vector<u_int8_t> *buff, const char *output_file,
                           map<u_int8_t, pair<bitset<MAX_COUNT_BITS_OF_CHAR>, u_int32_t>> &table) {
        map<u_int16_t, u_long> alphabetWithFrequency;
        create_alphabet_frequency(buff, alphabetWithFrequency);

        if (alphabetWithFrequency.size() == 1) {
            for (auto elem: alphabetWithFrequency) {
                pair<bitset<MAX_COUNT_BITS_OF_CHAR>, u_int8_t> coding_string;
                coding_string.first[0] = true;
                coding_string.second = 1;
                table.insert({elem.first, coding_string});
            }
        } else {
            create_table(alphabetWithFrequency, table);
        }

        vector<bitset<CHAR_BIT>> code_of_string;
        u_int32_t countOfWroteBits = create_compressing_string(buff, table, code_of_string);
        cout << table.size() * (1 + 1 + MAX_COUNT_BITS_OF_CHAR / 8) + 4 << "\n";

        write_to_bin(output_file, countOfWroteBits, code_of_string, table);
    }

    class EncodedSymbol {
        u_int8_t symbol{};
        string codes;

    public:
        EncodedSymbol() = default;

        explicit EncodedSymbol(u_int8_t symbol, pair<bitset<MAX_COUNT_BITS_OF_CHAR>, u_int32_t> tab) {
            this->symbol = symbol;
            for (u_int32_t i = 0; i < tab.second; ++i) {
                this->codes += to_string(tab.first[i]);
            }
        }

        u_int8_t getSymbol() const {
            return symbol;
        }

        const string &getCodes() const {
            return codes;
        }
    };

    bool compare_code_symbol(EncodedSymbol elem1, EncodedSymbol elem2) {
        return elem1.getCodes() < elem2.getCodes();
    }

    void print_encoded_symbols_in_lexicoraphical_order(
            const map<u_int8_t, pair<bitset<MAX_COUNT_BITS_OF_CHAR>, u_int32_t>> &table) {
        typedef set<EncodedSymbol, decltype(&compare_code_symbol)> lexicographical_order;
        lexicographical_order lexicographical_order_table(compare_code_symbol);

        for (auto element:table) {
            auto tab = new EncodedSymbol(element.first, element.second);
            lexicographical_order_table.insert(*tab);
            delete tab;
        }

        for (const auto &elem: lexicographical_order_table) {
            cout << elem.getCodes() << " " << static_cast<u_int32_t>( elem.getSymbol()) << "\n";
        }
    }
}

void compression(const char *input, const char *output_file, bool needPrint) {
    map<u_int8_t, std::pair<std::bitset<MAX_COUNT_BITS_OF_CHAR>, u_int32_t>> table;
    ifstream input_stream(input, ios::binary | ios::in);
    char ch;
    vector<u_int8_t> buff;
    while (input_stream.get(ch)) {
        buff.push_back(static_cast<u_int8_t>(ch));
    }
    input_stream.close();
    if (buff.empty()) {
        cout << "0\n0\n0\n";
        ofstream out;
        out.open(output_file, ios::binary | ios::out);
        out.close();
        return;
    }
    create_dictionary(&buff, output_file, table);

    if (needPrint && !table.empty()) {
        print_encoded_symbols_in_lexicoraphical_order(table);
    }
}

namespace {
    void read_from_stream_compressing_file(
            const char *input_file_name,
            vector<bitset<CHAR_BIT>> &bits,
            vector<u_int8_t> &bits_int,
            u_int32_t &COUNT_BYTES_IN_FILE) {
        ifstream input_stream(input_file_name, ios::binary);
        char input;
        while (input_stream.get(input)) {
            bits.emplace_back(bitset<sizeof(input) * CHAR_BIT>(input));
            for (int8_t k = 7; k >= 0; k--) {
                bits_int.push_back(bits[bits.size() - 1][k]);
            }
        }
        input_stream.close();
        COUNT_BYTES_IN_FILE = bits.size();
    }

//первые 32 бита
    u_int32_t find_count_of_bits_in_compressing_string(const vector<bitset<CHAR_BIT>> &all_bits) {
        bitset<32> bits;
        for (size_t k = 0; k < BYTE_IN_COUNT_OF_CODING_STRING; ++k) {
            for (size_t i = 0; i < BITS_COUNT; ++i) {
                bits[BITS_COUNT * k + i] = all_bits[k][i];
            }
        }
        return bits.to_ulong();
    }

    void print_count_bytes_compressing_string(
            u_int32_t count_bits_string,
            u_int32_t &BYTES_FIRST_STRING,
            u_int32_t &BITS_FIRST_STRING,
            u_int32_t &COUNT_OF_UNIQUE_ELEMENT,
            const u_int32_t COUNT_BYTES_IN_FILE) {
        BYTES_FIRST_STRING = ceil(static_cast<double> (count_bits_string) / 8) + BYTE_IN_COUNT_OF_CODING_STRING;
        BITS_FIRST_STRING = BYTES_FIRST_STRING * 8;
        COUNT_OF_UNIQUE_ELEMENT = (COUNT_BYTES_IN_FILE - BYTES_FIRST_STRING) / BYTES_IN_ONE_LINE_OF_TABLE;
        cout << static_cast<u_int32_t>(ceil(static_cast<double> (count_bits_string) / 8))
             << "\n"; // размер закодированного
    }

    void read_compressiong_string(
            const vector<u_int8_t> &all_bits,
            const u_int32_t count_of_bits,
            vector<u_int8_t> &coding_string) {
        for (size_t i = BITS_IN_COUNT_OF_CODING_STRING;
             i < BITS_IN_COUNT_OF_CODING_STRING + count_of_bits; ++i) {
            coding_string.push_back(all_bits[i]);
        }
    }

    char read_the_char(const vector<bitset<CHAR_BIT>> &all_bits,
                       const u_int32_t shift,
                       const u_int32_t BYTES_FIRST_STRING) {
        return static_cast<char>(all_bits[BYTES_FIRST_STRING + shift * BYTES_IN_ONE_LINE_OF_TABLE].to_ulong());
    }

    u_long read_count_of_code_char(const vector<bitset<CHAR_BIT>> &all_bits,
                                   const u_int32_t shift,
                                   const u_int32_t BYTES_FIRST_STRING) {
        bitset<8> bits;
        for (size_t i = 0; i < 8; ++i) {
            bits[i] = all_bits[BYTES_FIRST_STRING + shift * BYTES_IN_ONE_LINE_OF_TABLE + 1][i];

        }
        return bits.to_ulong();
    }

    void read_encoded_symbol_table_bin(
            const vector<bitset<CHAR_BIT>> &bits,
            const vector<u_int8_t> &bits_int,
            map<u_int8_t, pair<bitset<MAX_COUNT_BITS_OF_CHAR>, u_int32_t>> &table,
            map<u_int8_t, vector<char>> &tableForHuff,
            const u_int32_t BYTES_FIRST_STRING,
            const u_int32_t BITS_FIRST_STRING,
            const u_int32_t COUNT_OF_UNIQUE_ELEMENT) {
        u_int32_t i = 0;

        while (i < COUNT_OF_UNIQUE_ELEMENT) {
            unsigned char symbol = read_the_char(bits, i, BYTES_FIRST_STRING);
            u_long count_of_code_char = read_count_of_code_char(bits, i, BYTES_FIRST_STRING);

            vector<char> coding_char;
            bitset<MAX_COUNT_BITS_OF_CHAR> bit;
            read_code_char(bits_int, count_of_code_char, &coding_char, bit, i, BITS_FIRST_STRING);
            table[symbol] = {bit, count_of_code_char};
            tableForHuff[symbol] = coding_char;

            ++i;
        }
    }

    void decompress_compressing_string(
            const map<u_int8_t, vector<char>> &tableForHuff,
            const u_int32_t count_bits_string,
            const vector<u_int8_t> &coding_string,
            vector<u_int8_t> &string) {
        auto btHuffman = new BTHuffman();

        for (auto element:tableForHuff) {
            btHuffman->insert_to_Huffman_by_path(element.first, &element.second);
        }

        size_t i = 0;
        while (i != count_bits_string) {
            i = btHuffman->print_symbol(i, &coding_string, &string);
        }
        delete btHuffman->getRoot();
        delete btHuffman;
    }

    void write_to_file_decompress_string(const char *output_file, const vector<u_int8_t> &string) {
        ofstream out;
        out.open(output_file, ios::binary | ios::out);
        for (u_int8_t elem: string) {
            out.write(reinterpret_cast<const char *>(&elem), 1);
        }
        out.close();
    }

    void write_to_file_if_compress_empty(
            const vector<u_int8_t> &bits_int,
            const char *output_file) {
        if (bits_int.empty()) {
            cout << "0\n0\n0";
            ofstream out;
            out.open(output_file, ios::binary | ios::out);
            out.close();
            return;
        }
    }

    void print_size_of_encoded_symbol(const u_int32_t COUNT_OF_UNIQUE_ELEMENT) {
        cout << COUNT_OF_UNIQUE_ELEMENT * BYTES_IN_ONE_LINE_OF_TABLE + BYTE_IN_COUNT_OF_CODING_STRING << "\n";
    }
}

void decompression(
        const char *input_file_name,
        const char *output_file,
        bool needPrint) {
    map<u_int8_t, std::pair<std::bitset<MAX_COUNT_BITS_OF_CHAR>, u_int32_t>> table;
    vector<bitset<CHAR_BIT>> bits;
    vector<u_int8_t> bits_int;

    u_int32_t COUNT_BYTES_IN_FILE = 0;
    read_from_stream_compressing_file(input_file_name, bits, bits_int, COUNT_BYTES_IN_FILE);
    write_to_file_if_compress_empty(bits_int, output_file);
    if (bits_int.empty()) {
        return;
    }
    u_int32_t BYTES_FIRST_STRING = 0;
    u_int32_t BITS_FIRST_STRING = 0;
    u_int32_t COUNT_OF_UNIQUE_ELEMENT = 0;
    u_int32_t count_bits_string = find_count_of_bits_in_compressing_string(bits);
    print_count_bytes_compressing_string(count_bits_string, BYTES_FIRST_STRING, BITS_FIRST_STRING,
                                         COUNT_OF_UNIQUE_ELEMENT, COUNT_BYTES_IN_FILE);

    vector<u_int8_t> coding_string;
    read_compressiong_string(bits_int, count_bits_string, coding_string);

    map<u_int8_t, vector<char>> tableForHuff;
    read_encoded_symbol_table_bin(bits, bits_int, table, tableForHuff, BYTES_FIRST_STRING, BITS_FIRST_STRING,
                                  COUNT_OF_UNIQUE_ELEMENT);

    vector<u_int8_t> string;
    decompress_compressing_string(tableForHuff, count_bits_string, coding_string, string);
    cout << string.size() << "\n";

    print_size_of_encoded_symbol(COUNT_OF_UNIQUE_ELEMENT);
    write_to_file_decompress_string(output_file, string);

    if (needPrint && !table.empty()) {
        print_encoded_symbols_in_lexicoraphical_order(table);
    }
}
