/*
    Project Pegasus
    Version 0.1.0

    Copyright 2024-2026, Christopher A. Lupp

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published
    by the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/

#include "pegasus/philote/vlm_discipline.h"

#include <iostream>
#include <memory>
#include <string>

#include <grpcpp/grpcpp.h>

int main(int argc, char* argv[]) {
    std::string address = "localhost:50051";
    if (argc > 1)
        address = argv[1];

    auto discipline = std::make_shared<pegasus::VlmDiscipline>();

    grpc::ServerBuilder builder;
    builder.AddListeningPort(address, grpc::InsecureServerCredentials());
    discipline->RegisterServices(builder);

    std::unique_ptr<grpc::Server> server(builder.BuildAndStart());
    std::cout << "VLM Philote server listening on " << address << std::endl;
    server->Wait();

    return 0;
}
