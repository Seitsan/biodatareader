"""
Анализ SAM-файла с поддержкой командной строки. (Argparse-версия)

Example:
    # Получить заголовки
    python run_sam_argparse.py sample.sam header

    # Посчитать выравнивания
    python run_sam_argparse.py sample.sam count

    # Получить статистику
    python run_sam_argparse.py sample.sam stats
    
    # Запуск с фильтрацией по региону
    python run_sam_argparse.py sample.sam filter chr1 10000 20000
"""

import sys
from pathlib import Path
import argparse
from sam_reader import SamReader  # Предполагается, что этот модуль у вас есть


def handle_header(reader: SamReader):
    """1. Заголовки"""
    print("\n=== Заголовки SAM-файла ===")
    header = reader.get_header()
    if header:
        for tag, entries in header.items():
            print(f"{tag}:")
            # Убираем дубликаты и сортируем
            for line in sorted(set(entries)):
                print(f"  {line}")
    else:
        print("Заголовки отсутствуют.")


def handle_count(reader: SamReader):
    """2. Количество выравниваний"""
    total = reader.count_alignments()
    print(f"\n=== Общее количество выравниваний: {total:,}")


def handle_stats(reader: SamReader):
    """3. Статистика по хромосомам"""
    print("\n=== Статистика по хромосомам ===")
    df_stats = reader.stats_by_chromosome()
    if df_stats.empty:
        print("Нет выравниваний для анализа.")
    else:
        print(df_stats.to_string(index=False))


def handle_filter(reader: SamReader, chrom: str, start: int, end: int):
    """4. Фильтрация по региону"""
    if start > end:
        print("Ошибка: START не может быть больше END.", file=sys.stderr)
        sys.exit(1)

    print(f"\n=== Выравнивания в регионе {chrom}:{start}-{end} ===")
    found = False
    for rec in reader.filter_by_region(chrom, start, end):
        # Предполагаем, что rec имеет атрибуты id, chrom, start, end, cigar
        print(f"{rec.id}\t{rec.chrom}\t{rec.start}\t{rec.end}\t{rec.cigar}")
        found = True
    if not found:
        print("Нет выравниваний в указанном регионе.")


def main() -> None:
    """
    Основная функция CLI-утилиты для анализа SAM-файлов.
    """
    parser = argparse.ArgumentParser(
        description="Анализ SAM-файла: заголовки, статистика, фильтрация по региону."
    )
    
    # Глобальный аргумент
    parser.add_argument(
        "sam_file",
        type=Path,
        help="Путь к SAM-файлу"
    )

    # Контейнер для подкоманд
    subparsers = parser.add_subparsers(
        title="Команды",
        dest="command",
        required=True,
        help="Действие для выполнения"
    )

    # Команда "header"
    subparsers.add_parser(
        "header",
        help="Показать заголовки SAM-файла"
    )

    # Команда "count"
    subparsers.add_parser(
        "count",
        help="Посчитать общее число выравниваний"
    )

    # Команда "stats"
    subparsers.add_parser(
        "stats",
        help="Собрать статистику по хромосомам"
    )

    # Команда "filter"
    filter_parser = subparsers.add_parser(
        "filter",
        help="Фильтровать выравнивания по геномному региону"
    )
    filter_parser.add_argument(
        "chrom",
        help="Хромосома для фильтрации (например: '1', 'chr1')"
    )
    filter_parser.add_argument(
        "start",
        type=int,
        help="Начало региона (1-based)"
    )
    filter_parser.add_argument(
        "end",
        type=int,
        help="Конец региона"
    )

    args = parser.parse_args()

    # Проверка существования файла
    if not args.sam_file.is_file():
        print(f"Ошибка: файл не найден — {args.sam_file}", file=sys.stderr)
        sys.exit(1)

    # === Основной анализ ===
    try:
        with SamReader(args.sam_file) as reader:
            
            # Диспетчеризация команд
            if args.command == "header":
                handle_header(reader)
            
            elif args.command == "count":
                handle_count(reader)

            elif args.command == "stats":
                handle_stats(reader)

            elif args.command == "filter":
                handle_filter(reader, args.chrom, args.start, args.end)

    except Exception as e:
        print(f"Ошибка при обработке файла: {e}", file=sys.stderr)
        sys.exit(1)

    print()


if __name__ == "__main__":
    main()
