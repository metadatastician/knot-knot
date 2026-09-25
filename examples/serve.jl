# SPDX-License-Identifier: MPL-2.0
# serve.jl — stdlib-only static server for the web laboratory.
#
# Estate language policy keeps this free of Python/Node: a small blocking
# HTTP/1.1 file server over Sockets, rooted at www/public.
#
# Usage: julia --project=. examples/serve.jl [port]   (default 8081)

using Sockets

const ROOT = normpath(joinpath(@__DIR__, "..", "www", "public"))
const PORT = isempty(ARGS) ? 8081 : parse(Int, ARGS[1])

const MIME = Dict(
    ".html" => "text/html; charset=utf-8",
    ".css" => "text/css; charset=utf-8",
    ".js" => "text/javascript; charset=utf-8",
    ".json" => "application/json; charset=utf-8",
    ".svg" => "image/svg+xml",
    ".png" => "image/png",
    ".ico" => "image/x-icon",
    ".txt" => "text/plain; charset=utf-8",
)

function respond(client, status, body, ctype)
    write(client, "HTTP/1.1 $status\r\n")
    write(client, "Content-Type: $ctype\r\n")
    write(client, "Content-Length: $(sizeof(body))\r\n")
    write(client, "Connection: close\r\n\r\n")
    write(client, body)
end

function handle(client)
    line = readline(client)
    while !eof(client)
        l = readline(client)
        isempty(strip(l)) && break
    end
    parts = split(line, " ")
    path = length(parts) >= 2 ? String(parts[2]) : "/"
    path = split(path, "?")[1]
    rel = replace(path, ".." => "")
    rel = lstrip(rel, '/')
    rel = isempty(rel) ? "knot-lab/index.html" : rel
    file = normpath(joinpath(ROOT, rel))
    if !startswith(file, ROOT) || !isfile(file)
        respond(client, "404 Not Found", "not found", "text/plain; charset=utf-8")
    else
        ext = lowercase(splitext(file)[2])
        respond(client, "200 OK", read(file), get(MIME, ext, "application/octet-stream"))
    end
    close(client)
end

server = listen(ip"0.0.0.0", PORT)
println("knot-knot laboratory: http://0.0.0.0:$PORT  (root: $ROOT)")
while true
    client = accept(server)
    @async try
        handle(client)
    catch err
        @warn "serve: handler error" err
    end
end
