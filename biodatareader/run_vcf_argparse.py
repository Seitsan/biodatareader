"""
Демонстрационная программа для анализа VCF-файлов.

1. Получение заголовка и информации по отдельным группам заголовков (##INFO, ##FILTER и т.д.)
2. Получение количества вариантов.
3. Получение статистики “количество вариантов — регион” (регион = хромосома) с использованием pandas.
4. Получение вариантов в заданном геномном отрезке (аналог bedtools intersect).


Example:
    Запуск без региона:

    .. code-block:: bash

        python run_vcf_argparse.py sample.vcf

    Запуск с фильтрацией по региону:

    .. code-block:: bash

        python run_vcf_argparse.py sample.vcf chr1 10000 20000

"""

import sys
import argparse
from pathlib import Path
from vcf_reader import VcfReader

def handle_header(reader: VcfReader):
    """Выводит заголовки и группы заголовков."""
    print("=" * 70)
    print("1. ЗАГОЛОВКИ VCF-ФАЙЛА")
    print("=" * 70)
    header = reader.get_header()
    if header:
        print(f"Найдено {len(header)} мета-заголовков (##...)")
        for line in header[:5]:
            print(f"  {line}")
        if len(header) > 5:
            print(f"  ... и ещё {len(header) - 5} строк")
    else:
        print("Мета-заголовки не найдены.")

    print("\n" + "=" * 70)
    print("2. ИНФОРМАЦИЯ ПО ГРУППАМ ЗАГОЛОВКОВ")
    print("=" * 70)
    groups_to_show = ["INFO", "FILTER", "FORMAT", "contig"]
    for group in groups_to_show:
        entries = reader.get_header_group(group)
        if entries:
            print(f"\n##{group} — найдено записей: {len(entries)}")
            for entry in entries[:2]:  # Показываем максимум 2
                print(f"  {entry}")
            if len(entries) > 2:
                print(f"    ... и ещё {len(entries) - 2}")
        else:
            print(f"\n##{group} — не найдены")


def handle_count(reader: VcfReader):
    """Выводит общее количество вариантов."""
    print("=" * 70)
    print("3. КОЛИЧЕСТВО ВАРИАНТОВ")
    print("=" * 70)
    total = reader.count_variants()
    print(f"Общее количество вариантов: {total:,}")


def handle_stats(reader: VcfReader):
    """Выводит статистику по регионам (хромосомам)."""
    print("=" * 70)
    print("4. СТАТИСТИКА ПО РЕГИОНАМ (ХРОМОСОМАМ)")
    print("=" * 70)
    stats_df = reader.stats_by_region()
    if stats_df.empty:
        print("Нет вариантов для анализа.")
    else:
        print(stats_df.to_string(index=False))
        print(f"\nВсего регионов (хромосом) с вариантами: {len(stats_df)}")


def handle_filter(reader: VcfReader, chrom: str, start: int, end: int):
    """Выводит варианты в заданном регионе."""
    
    # Валидация региона
    if start > end or start < 1:
        print("Ошибка: START должен быть ≥ 1 и ≤ END.", file=sys.stderr)
        sys.exit(1)

    print(f"\n" + "=" * 70)
    print(f"5. ВАРИАНТЫ В РЕГИОНЕ: {chrom}:{start}-{end}")
    print("=" * 70)
    
    variants_in_region = list(reader.filter_by_region(chrom, start, end))
    print(f"Найдено вариантов: {len(variants_in_region)}")
    for i, var in enumerate(variants_in_region[:5], 1):  # Показываем первые 5
        print(f"  {i}. {var.chrom}:{var.pos} {var.ref}>{var.alt}")
    if len(variants_in_region) > 5:
        print(f"  ... и ещё {len(variants_in_region) - 5}")

def main() -> None:
    """
    Основная функция CLI-утилиты для анализа VCF-файлов.


    - Выводит мета-заголовки VCF-файла (строки, начинающиеся с '##').
    - Группирует и отображает информацию по ключевым секциям заголовка: ##INFO, ##FILTER, ##FORMAT, ##contig.
    - Подсчитывает общее количество вариантов в файле.
    - Формирует и выводит статистику по хромосомам (регионам) с использованием pandas.
    - (Опционально) Фильтрует и отображает варианты в заданном геномном регионе.

    Аргументы командной строки:

    - vcf_file (str): Обязательный путь к VCF-файлу.
    - chrom (str, optional): Название хромосомы для фильтрации (например, 'chr1').
    - start (int, optional): Начало региона (1-based, включительно).
    - end (int, optional): Конец региона (включительно).

    Программа завершается с кодом 1 в следующих случаях:

    - Указанный файл не существует.
    - Задана хромосома, но не указаны обе координаты (start и end).
    - Нарушены ограничения на координаты (start < 1 или start > end).
    - Произошла ошибка при чтении или обработке VCF-файла.

    Вывод:

    - Краткая сводка по заголовкам.
    - Статистика по группам метаинформации.
    - Общее число вариантов.
    - Таблица распределения вариантов по хромосомам.
    - Список вариантов в указанном регионе (максимум 5 первых, если их больше).

    Note:
        Для корректной работы требуется, чтобы класс VcfReader реализовывал
        следующие методы:
            - get_header() → list[str]
            - get_header_group(group: str) → list[str]
            - count_variants() → int
            - stats_by_region() → pandas.DataFrame
            - filter_by_region(chrom: str, start: int, end: int) → Iterator[VariantRecord]

    Raises:
        SystemExit: При ошибках валидации входных данных или обработки файла.
    """
    parser = argparse.ArgumentParser(
        description="Анализ VCF-файлов с использованием подкоманд.",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    # 2. Добавляем ГЛОБАЛЬНЫЙ аргумент (vcf_file), общий для всех
    parser.add_argument(
        "vcf_file",
        type=Path,
        help="Путь к входному VCF-файлу"
    )

    # 3. Добавляем контейнер для подкоманд
    subparsers = parser.add_subparsers(
        title="Команды",
        dest="command",  # Имя подкоманды будет сохранено в args.command
        required=True,   # Пользователь обязан выбрать одну из команд
        help="Действие для выполнения"
    )

    # 4. Подкоманда "header"
    subparsers.add_parser(
        "header",
        help="1. Показать заголовки (meta и группы INFO/FILTER/...)"
    )

    # 5. Подкоманда "count"
    subparsers.add_parser(
        "count",
        help="2. Посчитать общее количество вариантов"
    )

    # 6. Подкоманда "stats"
    subparsers.add_parser(
        "stats",
        help="3. Показать статистику по регионам (хромосомам)"
    )

    # 7. Подкоманда "filter" (со своими аргументами)
    filter_parser = subparsers.add_parser(
        "filter",
        help="4. Фильтровать варианты в заданном регионе"
    )
    filter_parser.add_argument(
        "chrom",
        help="Хромосома для фильтрации (например: '1', 'chr1')"
    )
    filter_parser.add_argument(
        "start",
        type=int,
        help="Начало региона (1-based, включительно)"
    )
    filter_parser.add_argument(
        "end",
        type=int,
        help="Конец региона (включительно)"
    )

    

    args = parser.parse_args()


    if not args.vcf_file.is_file():
        print(f"Ошибка: файл не найден — {args.vcf_file}", file=sys.stderr)
        sys.exit(1)

    try:
         with VcfReader(args.vcf_file) as reader:
            
            if args.command == "header":
                handle_header(reader)
            
            elif args.command == "count":
                handle_count(reader)
            
            elif args.command == "stats":
                handle_stats(reader)
            
            elif args.command == "filter":
                handle_filter(reader, args.chrom, args.start, args.end)

    except Exception as e:
        print(f"Ошибка при обработке VCF-файла: {e}", file=sys.stderr)
        sys.exit(1)

    print("\nАнализ завершён.")


if __name__ == "__main__":
    main()
