import pathlib
from contextlib import redirect_stdout
import click
import sys
import json
import csv
from fasta_reader import FastaReader

def output_results(data, output_format, command_name, output_file=None):
    if output_file:
        with open(output_file, 'w', encoding='utf-8') as f:
            with redirect_stdout(f):
                _output_data(data, output_format, command_name)
        print(f"Результаты сохранены в {output_file}")
    else:
        _output_data(data, output_format, command_name)


def _output_data(data, output_format, command_name):
    if output_format == "text":
        output_text(data, command_name)
    elif output_format == "json":
        output_json(data)
    elif output_format == "csv":
        output_csv(data, command_name)


def output_text(data, command_name):
    if command_name == "stats":
        print(f"Статистика FASTA-файла:")
        print(f"    Количество последовательностей: {data['count']}")
        print(f"    Средняя длина: {data['mean_length']:.2f}")

    elif command_name == "sequences":
        print(f"Найдено последовательностей: {len(data)}")
        for seq in data:
            print(f"    ID: {seq['id']}")
            print(f"    Длина: {seq['length']}")
            print("  " + "-" * 40)


def output_json(data):
    json.dump(data, sys.stdout, indent=2, ensure_ascii=False)
    print()


def output_csv(data, command_name):
    writer = csv.writer(sys.stdout)

    if command_name == "stats":
        writer.writerow(["metrics", "value"])
        writer.writerow(["Количество последовательностей", data['count']])
        writer.writerow(["Средняя длина", f"{data['mean_length']:.2f}"])

    elif command_name == "sequences":
        writer.writerow(["ID", "length"])
        for seq in data:
            writer.writerow([seq['id'], seq['length']])


@click.group()
def cli():
    pass

@cli.command()
@click.argument("filename", type=click.Path(exists=True, dir_okay=False, readable=True, path_type=pathlib.Path))
@click.option("-f", "--format", type=click.Choice(['text', 'json', 'csv']), default='text', help="Формат вывода")
def stats(filename:pathlib.Path, format):
    """Общая статистика"""
    try:
        with FastaReader(filename) as reader:
            list(reader.read())
            count = reader.get_seq_score()
            mean_len = reader.get_mean_seq_length()

            data = {
                'count': count,
                'mean_length': mean_len
            }

            output_results(data, format, "stats")
    except ValueError as e:
        print(f"Ошибка валидации последовательности: {e}", file=sys.stderr)
    except Exception as e:
        print(f"Ошибка при чтении файла: {e}", file=sys.stderr)


@cli.command()
@click.argument("filename", type=click.Path(exists=True, dir_okay=False, readable=True, path_type=pathlib.Path))
@click.option("-f", "--format", type=click.Choice(['text', 'json', 'csv']), default='text', help="Формат вывода")
@click.option("-o", "--output", type=click.Path(path_type=pathlib.Path), help="Путь к файлу для сохранения результатов")
def sequences(filename:pathlib.Path, format, output):
    """Работа с отдельными последовательностями"""

    sequences_data = []
    with FastaReader(filename) as reader:
        for record in reader.read():
            sequences_data.append({
                'id': record.id,
                'length': len(record.sequence),
                'sequence': record.sequence if format == "json" else None
            })

    output_results(sequences_data, format, "sequences", output)

if __name__ == "__main__":
    cli()