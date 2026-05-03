// Copyright (c) 2024 VLM Authors
// SPDX-License-Identifier: AGPL-3.0-or-later

#include <iostream>
#include <memory>
#include <string>

#include <grpcpp/grpcpp.h>

#include "vlm/philote/vlm_discipline.h"

int main(int argc, char* argv[]) {
  std::string address = "localhost:50051";
  if (argc > 1) address = argv[1];

  auto discipline = std::make_shared<vlm::VlmDiscipline>();

  grpc::ServerBuilder builder;
  builder.AddListeningPort(address, grpc::InsecureServerCredentials());
  discipline->RegisterServices(builder);

  std::unique_ptr<grpc::Server> server(builder.BuildAndStart());
  std::cout << "VLM Philote server listening on " << address << std::endl;
  server->Wait();

  return 0;
}
