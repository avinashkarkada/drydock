import click


@click.group()
@click.version_option()
def main() -> None:
    """Drydock: reproducible high-throughput virtual screening."""


if __name__ == "__main__":
    main()
