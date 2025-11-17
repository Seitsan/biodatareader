"""
Анализ SAM-файла с поддержкой командной строки. (Click-версия)

Example:
    # Получить заголовки
    python run_sam_click.py sample.sam header

    # Посчитать выравнивания
    python run_sam_click.py sample.sam count

    # Получить статистику
    python run_sam_click.py sample.sam stats
    
    # Запуск с фильтрацией по региону
    python run_sam_click.py sample.sam filter chr1 10000 20000
"""

import sys
from pathlib import Path
import click
from sam_reader import SamReader  # Предполагается, что этот модуль у вас есть


@click.group(context_settings=dict(help_option_names=['-h', '--help']))
@click.argument(
    "sam_file",
    type=click.Path(
        exists=True, 
        file_okay=True, 
        dir_okay=False, 
        readable=True, 
        path_type=Path
    )
)
@click.pass_context
def cli(ctx, sam_file: Path):
    """
    Анализ SAM-файла: заголовки, статистика, фильтрация по региону.
    
    SAM_FILE - Путь к входному SAM-файлу.
    """
    # Сохраняем путь в контексте для использования подкомандами
    ctx.obj = {"SAM_FILE": sam_file}


@cli.command()
@click.pass_context
def header(ctx):
    """1. Показать заголовки SAM-файла."""
    file_path = ctx.obj["SAM_FILE"]
    try:
        with SamReader(file_path) as reader:
            print("\n=== Заголовки SAM-файла ===")
            header_data = reader.get_header()
            if header_data:
                for tag, entries in header_data.items():
                    print(f"{tag}:")
                    # Убираем дубликаты и сортируем
                    for line in sorted(set(entries)):
                        print(f"  {line}")
            else:
                print("Заголовки отсутствуют.")
                
    except Exception as e:
        click.echo(f"Ошибка при обработке файла: {e}", err=True)
        sys.exit(1)
    print()


@cli.command()
@click.pass_context
def count(ctx):
    """2. Посчитать общее число выравниваний."""
    file_path = ctx.obj["SAM_FILE"]
    try:
        with SamReader(file_path) as reader:
            total = reader.count_alignments()
            print(f"\n=== Общее количество выравниваний: {total:,}")

    except Exception as e:
        click.echo(f"Ошибка при обработке файла: {e}", err=True)
        sys.exit(1)
    print()


@cli.command()
@click.pass_context
def stats(ctx):
    """3. Собрать статистику по хромосомам."""
    file_path = ctx.obj["SAM_FILE"]
    try:
        with SamReader(file_path) as reader:
            print("\n=== Статистика по хромосомам ===")
            df_stats = reader.stats_by_chromosome()
            if df_stats.empty:
                print("Нет выравниваний для анализа.")
            else:
                print(df_stats.to_string(index=False))

    except Exception as e:
        click.echo(f"Ошибка при обработке файла: {e}", err=True)
        sys.exit(1)
    print()


@cli.command()
@click.argument("chrom", type=str)
@click.argument("start", type=int)
@click.argument("end", type=int)
@click.pass_context
def filter(ctx, chrom: str, start: int, end: int):
    """4. Фильтровать выравнивания по геномному региону."""
    file_path = ctx.obj["SAM_FILE"]
    
    if start > end:
        click.echo("Ошибка: START не может быть больше END.", err=True)
        sys.exit(1)
        
    try:
        with SamReader(file_path) as reader:
            print(f"\n=== Выравнивания в регионе {chrom}:{start}-{end} ===")
            found = False
            for rec in reader.filter_by_region(chrom, start, end):
                # Предполагаем, что rec имеет атрибуты id, chrom, start, end, cigar
                print(f"{rec.id}\t{rec.chrom}\t{rec.start}\t{rec.end}\t{rec.cigar}")
                found = True
            if not found:
                print("Нет выравниваний в указанном регионе.")

    except Exception as e:
        click.echo(f"Ошибка при обработке файла: {e}", err=True)
        sys.exit(1)
    print()


if __name__ == "__main__":
    cli()
