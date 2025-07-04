{
  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs";
    flake-utils.url = "github:numtide/flake-utils";
  };

  outputs = { self, nixpkgs, flake-utils }:
    flake-utils.lib.eachDefaultSystem (system: let
      pkgs = import nixpkgs {
        inherit system;
      };
    in {
      devShells.default = pkgs.mkShell {
        buildInputs = with pkgs; [
          gcc
          cmake
          gdb

          clang-tools
          bear

          gnuplot
          jupyter

          (python3.withPackages(ps: with ps; [
            pybind11
          ]))
        ];
      };
    }
  );
}
