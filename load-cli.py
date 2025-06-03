import cli_funcs

netfabb_cone = cli_funcs.readCLI("/CAD/cli/cone_1.cli")
foss_file = cli_funcs.readCLI("/CAD/cli/foss-file.cli")
for layer in netfabb_cone[0:5]:
    print(layer.z)
for layer in foss_file[0:5]:
    print(layer.z)