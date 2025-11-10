"""
Анализ SAM-файла с поддержкой командной строки.

Example:
    Запуск без региона:

    .. code-block:: bash

        python run_sam_cl.py sample.sam

    Запуск с фильтрацией по региону:

    .. code-block:: bash

        python run_sam_cl.py sample.sam --chrom chr1 --start 10000 --end 20000
"""

from pathlib import Path
import click
from sam_reader import SamReader


@click.command()
@click.argument(
    "sam_file",
    type=click.Path(exists=True, path_type=Path, dir_okay=False),
    required=True
)
@click.option(
    "--chrom", "-c",
    type=str,
    help="Хромосома для фильтрации (например: '1', 'chr1')"
)
@click.option(
    "--start", "-s",
    type=int,
    help="Начало региона (1-based)"
)
@click.option(
    "--end", "-e",
    type=int,
    help="Конец региона"
)
@click.option(
    "--verbose", "-v",
    is_flag=True,
    help="Подробный вывод с дополнительной информацией"
)
def main(sam_file: Path, chrom: str, start: int, end: int, verbose: bool) -> None:
    """
    Анализ SAM-файла: заголовки, статистика, фильтрация по региону.

    SAM_FILE - Путь к SAM-файлу для анализа

    Примеры использования:
      python run_sam_cl.py sample.sam
      python run_sam_cl.py sample.sam --chrom chr1 --start 10000 --end 20000
      python run_sam_cl.py sample.sam -c chr1 -s 10000 -e 20000 -v
    """
    # Проверка корректности региона
    if chrom is not None:
        if start is None or end is None:
            raise click.BadParameter("Если задана хромосома, должны быть указаны и --start, и --end.")
        if start > end:
            raise click.BadParameter("START не может быть больше END.")

    # === Основной анализ ===
    try:
        with SamReader(sam_file) as reader:
            if verbose:
                click.echo("Начинаем анализ SAM-файла...")

            # === 1. Заголовки ===
            click.echo("\n" + click.style("=== Заголовки SAM-файла ===", bold=True))
            header = reader.get_header()
            if header:
                for tag, entries in header.items():
                    click.echo(click.style(f"{tag}:"))
                    for line in sorted(set(entries)):
                        click.echo(f"  {line}")
                if verbose:
                    click.echo(f"Найдено {len(header)} типов заголовков")
            else:
                click.echo("Заголовки отсутствуют.")

            # === 2. Количество выравниваний ===
            total = reader.count_alignments()
            click.echo("\n" + click.style(f"=== Общее количество выравниваний: {total:,} ===", bold=True))

            # === 3. Статистика по хромосомам ===
            click.echo("\n" + click.style("=== Статистика по хромосомам ===", bold=True))
            df_stats = reader.stats_by_chromosome()
            if df_stats.empty:
                click.echo("Нет выравниваний для анализа.")
            else:
                click.echo(df_stats.to_string(index=False))
                if verbose:
                    total_chroms = len(df_stats)
                    max_aligns = df_stats['count'].max()
                    min_aligns = df_stats['count'].min()
                    click.echo(f"Статистика: {total_chroms} хромосом, макс: {max_aligns:,}, мин: {min_aligns:,}")

            # === 4. Фильтрация по региону ===
            if chrom is not None:
                click.echo("\n" + click.style(f"=== Выравнивания в регионе {chrom}:{start}-{end} ===", fg='blue', bold=True))
                found_count = 0
                for rec in reader.filter_by_region(chrom, start, end):
                    click.echo(f"{rec.id}\t{rec.chrom}\t{rec.start}\t{rec.end}\t{rec.cigar}")
                    found_count += 1
                
                if found_count == 0:
                    click.echo("Нет выравниваний в указанном регионе.")
                elif verbose:
                    click.echo(f"Найдено выравниваний: {found_count}")

            if verbose:
                click.echo("\n" + click.style("Анализ завершен успешно!", bold=True))

    except Exception as e:
        click.echo(f"Ошибка при обработке файла: {e}", err=True)
        raise click.Abort()


if __name__ == "__main__":
    main()