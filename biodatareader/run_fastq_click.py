import click
from pathlib import Path
from analyze_fastq import analyze_fastq
from fastq_reader import FastqReader


@click.group()
def cli():
    """FASTQ Analyzer - инструмент для анализа FASTQ-файлов."""
    pass


def calculate_basic_stats(file_path):
    """Вычисляет базовую статистику по FASTQ-файлу."""
    total_sequences = 0
    total_length = 0
    sequence_lengths = []
    gc_count = 0
    total_bases = 0
    
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
                click.echo(f"Обработано записей: {total_sequences}", nl=False)
                click.echo("\r", nl=False)
    
    return {
        'total_sequences': total_sequences,
        'total_length': total_length,
        'sequence_lengths': sequence_lengths,
        'gc_count': gc_count,
        'total_bases': total_bases
    }


@cli.command()
@click.argument('file', type=click.Path(exists=True, path_type=Path))
@click.option('--detailed', '-d', is_flag=True, help='Подробная статистика с распределением длин')
def stats(file, detailed):
    """
    Вывод общей статистики по FASTQ-файлу.
    
    FILE - путь к FASTQ-файлу (поддерживается сжатие .gz)
    
    Выводит базовую статистику без построения графиков:
    - Общее количество последовательностей
    - Общее количество нуклеотидов  
    - Средняя/минимальная/максимальная длина
    - GC content
    """
    if not file.exists():
        raise click.ClickException(f"Файл не найден: {file}")
    
    click.echo(f"Статистика файла: {file}")
    click.echo("=" * 50)
    
    try:
        stats_data = calculate_basic_stats(file)
        
        if stats_data['total_sequences'] == 0:
            click.echo("В файле не найдено валидных последовательностей")
            return
        
        # Расчет статистики
        total_sequences = stats_data['total_sequences']
        total_length = stats_data['total_length']
        sequence_lengths = stats_data['sequence_lengths']
        gc_count = stats_data['gc_count']
        total_bases = stats_data['total_bases']
        
        mean_length = total_length / total_sequences
        min_length = min(sequence_lengths)
        max_length = max(sequence_lengths)
        gc_content = (gc_count / total_bases * 100) if total_bases > 0 else 0
        
        # Вывод основной статистики
        click.echo(f"{'Общее количество последовательностей:':<40} {total_sequences:>10,}")
        click.echo(f"{'Общее количество нуклеотидов:':<40} {total_length:>10,} bp")
        click.echo(f"{'Средняя длина последовательностей:':<40} {mean_length:>10.1f} bp")
        click.echo(f"{'Минимальная длина:':<40} {min_length:>10} bp")
        click.echo(f"{'Максимальная длина:':<40} {max_length:>10} bp")
        click.echo(f"{'GC content:':<40} {gc_content:>10.1f}%")
        
        # Детальная статистика по длинам
        if detailed:
            click.echo("\nРаспределение длин последовательностей:")
            length_counts = {}
            for length in sequence_lengths:
                length_counts[length] = length_counts.get(length, 0) + 1
            
            # Показываем топ-10 самых частых длин
            sorted_lengths = sorted(length_counts.items(), key=lambda x: x[1], reverse=True)[:10]
            for length, count in sorted_lengths:
                percentage = (count / total_sequences) * 100
                click.echo(f"  {length:>4} bp: {count:>8,} записей ({percentage:>5.1f}%)")
                
    except Exception as e:
        raise click.ClickException(f"Ошибка при анализе файла: {e}")


@cli.command()
@click.argument('file', type=click.Path(exists=True, path_type=Path))
def plots(file):
    """
    Построение графиков по FASTQ-файлу.
    
    FILE - путь к FASTQ-файлу (поддерживается сжатие .gz)
    
    Строит три графика:
    - Распределение длин последовательностей
    - Качество чтений по позициям
    - Содержание нуклеотидов по позициям
    """
    if not file.exists():
        raise click.ClickException(f"Файл не найден: {file}")
    
    click.echo(f"Построение графиков для файла: {file}")
    click.echo("-" * 50)
    
    try:
        # Сначала покажем базовую статистику
        stats_data = calculate_basic_stats(file)
        if stats_data['total_sequences'] == 0:
            click.echo("В файле не найдено валидных последовательностей")
            return
            
        click.echo(f"Обработано записей: {stats_data['total_sequences']:,}")
        click.echo("Строим графики...")
        
        # Запускаем построение графиков
        analyze_fastq(file)
        click.echo(" Графики успешно построены!")
        
    except Exception as e:
        raise click.ClickException(f"Ошибка при построении графиков: {e}")


@cli.command()
@click.argument('file', type=click.Path(exists=True, path_type=Path))
@click.option('--thresholds', '-t', is_flag=True, help='Показать пороги качества Q20/Q30')
def quality(file, thresholds):
    """
    Детальный анализ качества последовательностей.
    
    FILE - путь к FASTQ-файлу (поддерживается сжатие .gz)
    
    Анализирует качество чтений с расчетом:
    - Среднего/минимального/максимального качества
    - Процента оснований с Q≥20 и Q≥30
    - Общей оценки качества данных
    """
    if not file.exists():
        raise click.ClickException(f"Файл не найден: {file}")
    
    click.echo(f"Анализ качества для файла: {file}")
    click.echo("=" * 50)
    
    total_sequences = 0
    quality_stats = []
    
    try:
        with FastqReader(file) as reader:
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
                    click.echo(f"Обработано записей: {total_sequences}", nl=False)
                    click.echo("\r", nl=False)
        
        if total_sequences == 0:
            click.echo("В файле не найдено валидных последовательностей")
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
        
        # Вывод статистики качества
        click.echo(f"{'Общее количество последовательностей:':<40} {total_sequences:>10,}")
        click.echo(f"{'Общее количество нуклеотидов:':<40} {total_bases:>10,}")
        click.echo(f"{'Среднее качество (Phred):':<40} {mean_quality_all:>10.1f}")
        click.echo(f"{'Минимальное качество:':<40} {min_quality_all:>10}")
        click.echo(f"{'Максимальное качество:':<40} {max_quality_all:>10}")
        click.echo(f"{'Процент оснований с Q≥20:':<40} {q20_percentage:>9.1f}%")
        click.echo(f"{'Процент оснований с Q≥30:':<40} {q30_percentage:>9.1f}%")
        
        # Категоризация качества
        if q30_percentage >= 80:
            quality_grade = "Отличное"
            grade_color = "green"
        elif q30_percentage >= 60:
            quality_grade = "Хорошее" 
            grade_color = "blue"
        elif q20_percentage >= 80:
            quality_grade = "Удовлетворительное"
            grade_color = "yellow"
        else:
            quality_grade = "Плохое"
            grade_color = "red"
            
        click.echo(f"{'Оценка качества:':<40} ", nl=False)
        click.secho(f"{quality_grade:>10}", fg=grade_color, bold=True)
        
        # Дополнительная информация о порогах качества
        if thresholds:
            click.echo("\nПояснение по порогам качества:")
            click.echo("  Q20: 1% ошибок (99% точность) - минимальный приемлемый уровень")
            click.echo("  Q30: 0.1% ошибок (99.9% точность) - хорошее качество")
            click.echo("  Q40: 0.01% ошибок (99.99% точность) - отличное качество")
            
    except Exception as e:
        raise click.ClickException(f"Ошибка при анализе качества: {e}")


if __name__ == "__main__":
    cli()