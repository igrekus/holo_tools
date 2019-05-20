#include <iostream>
#include <string.h>

std::string read_image(std::string const &in_file_name) {
    std::cout << "read image..." << std::endl;
    return "image_data";
}

std::string encode_kinoform(std::string const &image) {
    std::cout << "apply Gerchberg-Saxton algorithm..." << std::endl;
    return "processed_data";
}

void write_image(std::string const &processed_data, std::string const &out_file_name) {
    std::cout << "write output file: " << out_file_name << std::endl;
}

void kin_prep(std::string const &in_file_name, std::string const &out_file_name) {
    std::cout << "preparing kinoform from image: " << in_file_name << std::endl;

    std::string image_data = read_image(in_file_name);
    std::string processed_data = encode_kinoform(image_data);
    write_image(processed_data, out_file_name);
}

void kin_view(std::string const &in_file_name) {
    std::cout << "reverting kinoform into an image: " << in_file_name << std::endl;
}

int main(int argc, char* argv[]) {

    std::string mode = argv[1];
    std::string in_file_name = argv[2];
    std::string out_file_name = argv[3];

    if (mode == "--kin-prep") {
        kin_prep(in_file_name, out_file_name);
    } else if (mode == "--kin-view") {
        kin_view(in_file_name);
    }

    return 0;
}