import argparse
import sys
from pathlib import Path
from analyze_fastq import analyze_fastq
from fastq_reader import FastqReader


def stats_command(args):
    """
    Выводит общую статистику по FASTQ-файлу без графиков.
    
    Args:
        args: Аргументы командной строки
    """
    file_path = Path(args.file)
    
    if not file_path.exists():
        print(f"Ошибка: файл не найден — {file_path}")
        sys.exit(1)
    
    print(f"Статистика файла: {file_path}")
    print("=" * 50)
    
    total_sequences = 0
    total_length = 0
    sequence_lengths = []
    gc_count = 0
    total_bases = 0
    
    try:
        with FastqReader(file_path) as reader:
            for record in reader.read():
                total_sequences += 1
                seq_len = len(record.sequence)
                total_length += seq_len
                sequence_lengths.append(seq_len)
                
                # Считаем GC content
                seq_upper = record.sequence.upper()
                gc_count += seq_upper.count('G') + seq_upper.count('C')
                total_bases += seq_len
                
                # Прогресс для больших файлов
                if total_sequences % 10000 == 0:
                    print(f"Обработано записей: {total_sequences}", end='\r')
        
        if total_sequences == 0:
            print("В файле не найдено валидных последовательностей")
            return
        
        # Расчет статистики
        mean_length = total_length / total_sequences
        min_length = min(sequence_lengths)
        max_length = max(sequence_lengths)
        gc_content = (gc_count / total_bases * 100) if total_bases > 0 else 0
        
        # Вывод статистики
        print(f"{'Общее количество последовательностей:':<40} {total_sequences:>10,}")
        print(f"{'Общее количество нуклеотидов:':<40} {total_length:>10,} bp")
        print(f"{'Средняя длина последовательностей:':<40} {mean_length:>10.1f} bp")
        print(f"{'Минимальная длина:':<40} {min_length:>10} bp")
        print(f"{'Максимальная длина:':<40} {max_length:>10} bp")
        print(f"{'GC content:':<40} {gc_content:>10.1f}%")
        
        # Дополнительная статистика по длинам
        print("\nРаспределение длин последовательностей:")
        length_counts = {}
        for length in sequence_lengths:
            length_counts[length] = length_counts.get(length, 0) + 1
        
        # Показываем топ-10 самых частых длин
        sorted_lengths = sorted(length_counts.items(), key=lambda x: x[1], reverse=True)[:10]
        for length, count in sorted_lengths:
            percentage = (count / total_sequences) * 100
            print(f"  {length} bp: {count:,} записей ({percentage:.1f}%)")
            
    except Exception as e:
        print(f"Ошибка при анализе файла: {e}")
        sys.exit(1)


def plots_command(args):
    """
    Строит графики по FASTQ-файлу.
    
    Args:
        args: Аргументы командной строки
    """
    file_path = Path(args.file)
    
    if not file_path.exists():
        print(f"Ошибка: файл не найден — {file_path}")
        sys.exit(1)
    
    print(f"Построение графиков для файла: {file_path}")
    print("-" * 50)
    
    try:
        analyze_fastq(file_path)
        print("Графики успешно построены!")
        
    except Exception as e:
        print(f"Ошибка при построении графиков: {e}")
        sys.exit(1)


def quality_command(args):
    """
    Детальный анализ качества последовательностей.
    
    Args:
        args: Аргументы командной строки
    """
    file_path = Path(args.file)
    
    if not file_path.exists():
        print(f"Ошибка: файл не найден — {file_path}")
        sys.exit(1)
    
    print(f"Анализ качества для файла: {file_path}")
    print("=" * 50)
    
    total_sequences = 0
    quality_stats = []
    
    try:
        with FastqReader(file_path) as reader:
            for record in reader.read():
                total_sequences += 1
                qual_scores = record.quality
                
                # Статистика по качеству для этой последовательности
                seq_quality = {
                    'mean_quality': sum(qual_scores) / len(qual_scores),
                    'min_quality': min(qual_scores),
                    'max_quality': max(qual_scores),
                    'q20_count': sum(1 for q in qual_scores if q >= 20),
                    'q30_count': sum(1 for q in qual_scores if q >= 30),
                    'length': len(qual_scores)
                }
                quality_stats.append(seq_quality)
                
                if total_sequences % 10000 == 0:
                    print(f"Обработано записей: {total_sequences}", end='\r')
        
        if total_sequences == 0:
            print("В файле не найдено валидных последовательностей")
            return
        
        # Агрегированная статистика по качеству
        total_bases = sum(stats['length'] for stats in quality_stats)
        total_q20 = sum(stats['q20_count'] for stats in quality_stats)
        total_q30 = sum(stats['q30_count'] for stats in quality_stats)
        
        mean_quality_all = sum(stats['mean_quality'] for stats in quality_stats) / total_sequences
        min_quality_all = min(stats['min_quality'] for stats in quality_stats)
        max_quality_all = max(stats['max_quality'] for stats in quality_stats)
        
        q20_percentage = (total_q20 / total_bases) * 100
        q30_percentage = (total_q30 / total_bases) * 100
        
        print(f"{'Общее количество последовательностей:':<40} {total_sequences:>10,}")
        print(f"{'Общее количество нуклеотидов:':<40} {total_bases:>10,}")
        print(f"{'Среднее качество (Phred):':<40} {mean_quality_all:>10.1f}")
        print(f"{'Минимальное качество:':<40} {min_quality_all:>10}")
        print(f"{'Максимальное качество:':<40} {max_quality_all:>10}")
        print(f"{'Процент оснований с Q≥20:':<40} {q20_percentage:>9.1f}%")
        print(f"{'Процент оснований с Q≥30:':<40} {q30_percentage:>9.1f}%")
        
        
    except Exception as e:
        print(f"Ошибка при анализе качества: {e}")
        sys.exit(1)

def main():
    """
    Основная функция с использованием argparse и подкоманд.
    """
    parser = argparse.ArgumentParser(
        description="FASTQ Analyzer - комплексный инструмент для анализа FASTQ-файлов",
        prog="run_fastq.py"
    )
    
    subparsers = parser.add_subparsers(
        title="доступные команды",
        description="выберите команду для анализа FASTQ-файлов",
        help="дополнительная информация по команде",
        dest="command",
        required=True
    )
    
    # Подкоманда stats - общая статистика
    stats_parser = subparsers.add_parser(
        "stats",
        help="вывод общей статистики по FASTQ-файлу",
        description="Вывод подробной статистики без построения графиков"
    )
    stats_parser.add_argument(
        "file",
        type=str,
        help="путь к FASTQ-файлу (поддерживается сжатие .gz)"
    )
    stats_parser.set_defaults(func=stats_command)
    
    # Подкоманда plots - построение графиков
    plots_parser = subparsers.add_parser(
        "plots",
        help="построение графиков по FASTQ-файлу",
        description="Визуализация данных через построение графиков"
    )
    plots_parser.add_argument(
        "file",
        type=str,
        help="путь к FASTQ-файлу (поддерживается сжатие .gz)"
    )
    plots_parser.set_defaults(func=plots_command)
    
    # Подкоманда quality - анализ качества
    quality_parser = subparsers.add_parser(
        "quality",
        help="детальный анализ качества последовательностей",
        description="Подробный анализ качества чтений с Q20/Q30 статистикой"
    )
    quality_parser.add_argument(
        "file",
        type=str,
        help="путь к FASTQ-файлу (поддерживается сжатие .gz)"
    )
    quality_parser.set_defaults(func=quality_command)
    

    # Парсинг аргументов
    args = parser.parse_args()
    
    # Вызов соответствующей функции
    args.func(args)


if __name__ == "__main__":
    main()