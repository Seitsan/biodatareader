"""
Демонстрационная программа для анализа VCF-файлов с использованием 
под-команд click.

Example:
    # Получить заголовки
    python run_vcf_click.py sample.vcf header

    # Посчитать варианты
    python run_vcf_click.py sample.vcf count

    # Получить статистику
    python run_vcf_click.py sample.vcf stats

    # Запустить фильтрацию по региону
    python run_vcf_click.py sample.vcf filter chr1 10000 20000
"""

import sys
import click
from pathlib import Path
from vcf_reader import VcfReader   

@click.group(context_settings=dict(help_option_names=['-h', '--help']))
@click.argument(
    "vcf_file",
    type=click.Path(
        exists=True,        
        file_okay=True,     
        dir_okay=False,     
        readable=True,      
        path_type=Path     
    ),
)
@click.pass_context  
def cli(ctx, vcf_file: Path):
    """
    Анализ VCF-файлов: заголовки, статистика, фильтрация.

    VCF_FILE - Путь к входному VCF-файлу.
    """
    ctx.obj = {"VCF_FILE": vcf_file}


@cli.command()
@click.pass_context   
def header(ctx):
    """1. Показать заголовки (meta и группы INFO/FILTER/...)."""
    
    file_path = ctx.obj["VCF_FILE"]  
    click.echo(f"--- Анализ заголовков для: {file_path.name} ---")

    try:
        with VcfReader(file_path) as reader:
            click.echo("=" * 70)
            click.echo("1. ОБЩИЕ ЗАГОЛОВКИ VCF-ФАЙЛА")
            click.echo("=" * 70)
            header_lines = reader.get_header()
            if header_lines:
                click.echo(f"Найдено {len(header_lines)} мета-заголовков (##...)")
                for line in header_lines[:5]:
                    click.echo(f"  {line}")
                if len(header_lines) > 5:
                    click.echo(f"  ... и ещё {len(header_lines) - 5} строк")
            else:
                click.echo("Мета-заголовки не найдены.")

            click.echo("\n" + "=" * 70)
            click.echo("2. ИНФОРМАЦИЯ ПО ГРУППАМ ЗАГОЛОВКОВ")
            click.echo("=" * 70)
            groups_to_show = ["INFO", "FILTER", "FORMAT", "contig"]
            for group in groups_to_show:
                entries = reader.get_header_group(group)
                if entries:
                    click.echo(f"\n##{group} — найдено записей: {len(entries)}")
                    for entry in entries[:2]:  # Показываем максимум 2
                        click.echo(f"  {entry}")
                    if len(entries) > 2:
                        click.echo(f"    ... и ещё {len(entries) - 2}")
                else:
                    click.echo(f"\n##{group} — не найдены")

    except Exception as e:
        click.echo(f"Ошибка при обработке VCF-файла: {e}", err=True)
        sys.exit(1)
    
    click.echo("\n--- Анализ заголовков завершён ---")


@cli.command()
@click.pass_context
def count(ctx):
    """2. Посчитать общее количество вариантов."""
    
    file_path = ctx.obj["VCF_FILE"]
    click.echo(f"--- Подсчет вариантов в: {file_path.name} ---")
    
    try:
        with VcfReader(file_path) as reader:
            click.echo("=" * 70)
            click.echo("3. КОЛИЧЕСТВО ВАРИАНТОВ")
            click.echo("=" * 70)
            total = reader.count_variants()
            click.echo(f"Общее количество вариантов: {total:,}")

    except Exception as e:
        click.echo(f"Ошибка при обработке VCF-файла: {e}", err=True)
        sys.exit(1)


@cli.command()
@click.pass_context
def stats(ctx):
    """3. Показать статистику по регионам (хромосомам)."""
    
    file_path = ctx.obj["VCF_FILE"]
    click.echo(f"--- Статистика по регионам в: {file_path.name} ---")

    try:
        with VcfReader(file_path) as reader:
            click.echo("=" * 70)
            click.echo("4. СТАТИСТИКА ПО РЕГИОНАМ (ХРОМОСОМАМ)")
            click.echo("=" * 70)
            stats_df = reader.stats_by_region()
            if stats_df.empty:
                click.echo("Нет вариантов для анализа.")
            else:
                click.echo(stats_df.to_string(index=False))
                click.echo(f"\nВсего регионов (хромосом) с вариантами: {len(stats_df)}")

    except Exception as e:
        click.echo(f"Ошибка при обработке VCF-файла: {e}", err=True)
        sys.exit(1)


@cli.command()
@click.argument("chrom", type=str)
@click.argument("start", type=int)
@click.argument("end", type=int)
@click.pass_context
def filter(ctx, chrom: str, start: int, end: int):
    """4. Фильтровать варианты в заданном регионе."""
    
    file_path = ctx.obj["VCF_FILE"]
    
     
    if start > end or start < 1:
        click.echo("Ошибка: START должен быть ≥ 1 и ≤ END.", err=True)
        sys.exit(1)

    click.echo(f"--- Фильтрация региона {chrom}:{start}-{end} в: {file_path.name} ---")

    try:
        with VcfReader(file_path) as reader:
            click.echo("=" * 70)
            click.echo(f"5. ВАРИАНТЫ В РЕГИОНЕ: {chrom}:{start}-{end}")
            click.echo("=" * 70)
            
            variants_in_region = list(reader.filter_by_region(chrom, start, end))
            
            click.echo(f"Найдено вариантов: {len(variants_in_region)}")
            for i, var in enumerate(variants_in_region[:5], 1):  
                 
                click.echo(f"  {i}. {var.chrom}:{var.pos} {var.ref}>{var.alt}")
            if len(variants_in_region) > 5:
                click.echo(f"  ... и ещё {len(variants_in_region) - 5}")

    except Exception as e:
        click.echo(f"Ошибка при обработке VCF-файла: {e}", err=True)
        sys.exit(1)


if __name__ == "__main__":
    cli()
