using Weave

jmd_files=filter(endswith("Jmd"),readdir(@__DIR__))

weave("$(@__DIR__)/Example.Jmd",doctype = "github", out_path = @__DIR__)
