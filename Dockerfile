# Drydock container image.
#
# Built from pixi.lock, exactly as a local install is, so the container and a
# workstation install cannot drift apart. CI runs the self-test in both, and a
# disagreement between them is a packaging bug.
#
# For HPC, convert to Apptainer rather than running Docker:
#
#     docker build -t drydock:0.1.0 .
#     apptainer build drydock.sif docker-daemon://drydock:0.1.0
#
# Apptainer, not Docker, is the right target for shared clusters: it is rootless
# and needs no daemon, whereas most sites refuse Docker precisely because access
# to the daemon is equivalent to root on the node.

FROM ghcr.io/prefix-dev/pixi:0.48.0 AS build

WORKDIR /app

# Copy only what the solve depends on first, so editing source does not
# invalidate the (slow) environment layer.
COPY pixi.toml pixi.lock pyproject.toml README.md ./
COPY src/ ./src/

# --locked aborts if pixi.lock has drifted from pixi.toml, so the build fails
# here rather than producing an image that quietly differs from what developers
# and CI are running. (It cannot be combined with --frozen, which only declines
# to update the lockfile rather than checking it.)
RUN pixi install --locked -e default

# `pixi shell-hook` bakes the activation into a plain script, so the runtime
# stage needs no pixi and no shell initialisation.
RUN pixi shell-hook -e default > /shell-hook.sh \
    && echo 'exec "$@"' >> /shell-hook.sh


FROM ubuntu:24.04 AS runtime

LABEL org.opencontainers.image.title="Drydock" \
      org.opencontainers.image.description="Reproducible high-throughput virtual screening" \
      org.opencontainers.image.licenses="GPL-3.0-or-later"

WORKDIR /app

COPY --from=build /app/.pixi/envs/default /app/.pixi/envs/default
COPY --from=build /shell-hook.sh /shell-hook.sh
COPY --from=build /app/src /app/src
COPY pyproject.toml README.md ./

# Runs are written to a directory the caller mounts; nothing of consequence
# lives inside the image.
VOLUME ["/data"]
WORKDIR /data

ENTRYPOINT ["/bin/bash", "/shell-hook.sh"]
CMD ["drydock", "--help"]
