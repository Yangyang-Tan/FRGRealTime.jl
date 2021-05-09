using Weave

jmd_files=filter(endswith("jmd"),readdir(@__DIR__))

jmd_path=joinpath.(@__DIR__,jmd_files)
html_path=joinpath(@__DIR__,"html")
markdown_path=joinpath(@__DIR__,"markdown")

weave.(jmd_path,doctype = "github", out_path = markdown_path)

weave.(jmd_path, out_path = html_path)
