import sys
import argparse
import json
import csv
from contextlib import redirect_stdout
from itertools import count
import pathlib
from fasta_reader import FastaReader


def setup_parser():
    parser = argparse.ArgumentParser(description="Инструмент для для анализа FASTA-файлов")

    parser.add_argument("-v", "--verbose", action="store_true", help="Подробный вывод")

    subparser = parser.add_subparsers(dest="command", help="Доступные команды")

    stats_parser = subparser.add_parser("stats", help="Общая статистика")
    stats_parser.add_argument("filename", type=pathlib.Path, help="Путь к FASTA-файлу")
    stats_parser.add_argument("--format", choices=["text", "json", "csv"], default="text", help="Формат вывода")

    sequences_parser = subparser.add_parser("sequences", help="Работа с отдельными последовательностями")
    sequences_parser.add_argument("filename", type=pathlib.Path, help="Путь к FASTA-файлу")
    sequences_parser.add_argument("--output", type=pathlib.Path, help="Путь к файлу для сохранения результатов")
    sequences_parser.add_argument("--format", choices=["text", "json", "csv"], default="text", help="Формат вывода")

    return parser


def output_results(data, output_format, command_name, output_file=None):

    if output_file:
        with open(output_file, 'w') as f:
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
        print(f"    Количество последовательностей: {len(data)}")
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
        writer.writerow(["Метрика", "Значение"])
        writer.writerow(["Количество последовательностей", data['count']])
        writer.writerow(["Средняя длина", f"{data['mean_length']:.2f}"])

    elif command_name == "sequences":
        writer.writerow(["ID", "Длина"])
        for seq in data:
            writer.writerow([seq['id'], seq['length']])


def handle_stats(args) -> None:

    try:
        with FastaReader(args.filename) as reader:

            list(reader.read())

            count = reader.get_seq_score()
            mean_len = reader.get_mean_seq_length()

            data = {
                'count': count,
                'mean_length': mean_len
            }

            output_results(data, args.format, "stats")

    except ValueError as e:
        print(f"Ошибка валидации последовательности: {e}", file=sys.stderr)
    except Exception as e:
        print(f"Ошибка при чтении файла: {e}", file=sys.stderr)


def handle_sequences(args):
    sequences_data = []
    with FastaReader(args.filename) as reader:
        for record in reader.read():
            sequences_data.append({
                'id': record.id,
                'length': len(record.sequence),
                'sequence': record.sequence if args.format == "json" else None
            })

    output_results(sequences_data, args.format, "sequences", args.output)

def main():
    parser = setup_parser()
    args = parser.parse_args()

    if args.command == "stats":
        handle_stats(args)
    elif args.command == "sequences":
        handle_sequences(args)


if __name__ == "__main__":
    main()