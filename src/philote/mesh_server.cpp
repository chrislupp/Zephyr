// Copyright (c) 2024 VLM Authors
// SPDX-License-Identifier: AGPL-3.0-or-later

#include "vlm/philote/mesh_discipline.h"

#include <iostream>
#include <memory>
#include <string>

#include <grpcpp/grpcpp.h>

int main(int argc, char* argv[]) {
    std::string address = "localhost:50052";
    if (argc > 1)
        address = argv[1];

    auto discipline = std::make_shared<vlm::MeshDiscipline>();

    grpc::ServerBuilder builder;
    builder.AddListeningPort(address, grpc::InsecureServerCredentials());
    discipline->RegisterServices(builder);

    std::unique_ptr<grpc::Server> server(builder.BuildAndStart());
    std::cout << "Mesh Philote server listening on " << address << std::endl;
    server->Wait();

    return 0;
}
