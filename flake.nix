{
  inputs = {
    nixpkgs.url = "github:nixos/nixpkgs?ref=nixos-unstable";
    devenv.url = "github:cachix/devenv";
    personal = { url = "sourcehut:~showyourcode/flakes"; inputs.nixpkgs.follows = "nixpkgs"; };
  };

  outputs = inputs@{ flake-parts, nixpkgs, ... }:
    flake-parts.lib.mkFlake { inherit inputs; } {
      imports = [
        inputs.devenv.flakeModule
      ];
      systems = nixpkgs.lib.systems.flakeExposed;

      perSystem = { config, self', inputs', pkgs, system, ... }: {
        # Per-system attributes can be defined here. The self' and inputs'
        # module parameters provide easy access to attributes of the same
        # system.
        _module.args.pkgs = import nixpkgs {
          inherit system;
          overlays = [
            inputs.personal.overlays.default
          ];
        };

        devenv.shells.default = {
          packages = with pkgs; [
            easel
            infernal
            nextflow
            graphviz
          ];

          languages = {
            python = {
              enable = true;
              poetry.enable = true;
            };

            typescript = {
              enable = true;
            };

            javascript = {
              enable = true;
              pnpm.enable = true;
            };
          };

          pre-commit.hooks = {
            black.enable = true;
            isort.enable = true;
            prettier.enable = true;
          };
        };
      };
    };
}
