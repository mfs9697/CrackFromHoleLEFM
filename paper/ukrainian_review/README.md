# IAM — авторське змістовне редагування українського рукопису

**Режим AAM-UA8 / IAM-UA: робоча авторська редакція, не дослівний переклад.** База — PR #91, коміт `57992bac6c747be71c4aa413af436e4137eb1869`. Повний текст, рисунки, бібліографія та аудити Claude збережені в базовій гілці. У цій гілці редагуємо погоджені розділи послідовно. Журнал рішень: [IAM_EDITORIAL_LEDGER.md](audit/IAM_EDITORIAL_LEDGER.md).

**Увага:** успадкований PDF та звіти про повне збереження перекладу і складання належать PR #91. Після змістовного редагування вони не підтверджують відповідність новому TeX. Перевірки дослівного перекладу не є критерієм прийнятності IAM-редакції; не послаблювати їх і не перезаписувати історичні звіти. Повторне складання та перевірка посилань потрібні для кожної погодженої редакції. Перевірку рисунків виконуємо разом із відповідними розділами.

---

# Українська редакція рукопису AAM

Робочий переклад для спільного наукового обговорення з професором Богдановим.
Джерело — `paper/submissions/aam/manuscript.tex` у коміті
`7597cba1a378e6f436a3bfab304f85702902ee4e` репозиторію CrackFromHoleLEFM.
Відповідність Git blob, SHA-256, історію PR і відмінності пізнішої гілки
зафіксовано в [аудиті джерела](audit/SOURCE_REVIEW.md).

[Успадкований PDF PR #91 (не відповідає новому TeX)](AAM_UA_Bogdanov_review.pdf) ·
[Початкові файли та ідентифікатори](SOURCE_MANIFEST.json) ·
[Наукові питання](audit/SCIENTIFIC_QUESTIONS.md) ·
[Звіт про рисунки](audit/FIGURE_REBUILD_REPORT.md) ·
[Мовний огляд](audit/CLAUDE_LINGUISTIC_REVIEW.md)

## Відтворення PDF (Windows, PowerShell)

Потрібні:

- XeLaTeX і BibTeX;
- STIX Two Text (Regular, Bold, Italic, Bold Italic) і STIX Two Math як файли `STIXTwoText-*.otf` і
  `STIXTwoMath-Regular.otf`, доступні XeLaTeX (пакет `stix2-otf`);
- пакети unicode-math, polyglossia, pgfplots і українські шаблони переносу.

Перевірено з MiKTeX 26.5. З кореня репозиторію:

```powershell
.\paper\ukrainian_review\build_ua.ps1
```

Сценарій виконує XeLaTeX → BibTeX → XeLaTeX → XeLaTeX. Перед копіюванням PDF він строго перевіряє журнал:
нерозв’язані посилання, переповнення, пропущені гліфи, шаблони переносу. Для збірки рукопису достатньо файлів
`main_ua.tex`, `references.bib` і `figures/`. Для сценаріїв аудиту потрібна копія репозиторію з історією, що
містить базовий коміт.

## Рисунки (AAM-UA7)

Усі 31 графічний файл у `figures/` створено командою:

```powershell
.\paper\ukrainian_review\scripts\build_ua_figures.ps1
```

Сценарій компілює панелі `figures_src/panels/**/*.tex` з XeLaTeX + PGFPlots і центральним стилем
`figures_src/ua_figstyle.tex`. Кожну панель спроєктовано за її остаточною друкованою шириною. Склад:

- 22 панелі побудовано з таблиць `figures_src/data/*.csv`;
- 9 панелей (рис. 1, 2(a–c), 8(a–d), 9(a)) містять незмінні векторні шари `figures_src/retained/`, на які наново
  набрано осі й написи.

Допоміжні сценарії потрібні лише для повного відтворення проміжних файлів. Вони не перезаписують наявних таблиць:

```powershell
# таблиці даних із записів (MATLAB R2023a, лише читання .mat), з кореня репозиторію
matlab -batch "run('paper/ukrainian_review/scripts/extract_ua_figure_data.m')"
matlab -batch "run('paper/ukrainian_review/scripts/dump_ua_mat_records.m')"
python paper/ukrainian_review/scripts/make_ua_constants.py
# збережені векторні шари з затверджених англійських PDF
python paper/ukrainian_review/scripts/build_retained_layers.py
```

## Історичні перевірки AAM-UA7 (для бази PR #91)

Не запускати команди з `--write-audit` на IAM-редакції до підготовки окремого редакторського протоколу: ці звіти належать повному перекладу.


```powershell
# незалежна перевірка даних рисунків і 9 негативних перевірок -> audit/FIGURE_DATA_VERIFICATION.json
python paper/ukrainian_review/scripts/run_ua_figure_verification.py
python paper/ukrainian_review/scripts/make_ua_figure_manifest.py
# збереження наукового змісту відносно закріплених Git-об'єктів і перевірка зібраного PDF
python paper/ukrainian_review/scripts/verify_preservation.py --write-audit
python paper/ukrainian_review/scripts/verify_review_pdf.py paper/ukrainian_review/.build --write-audit
```

Потрібен Python 3 з PyMuPDF і pypdf. Скінченно-елементних розрахунків немає. Збірку й перевірки виконано у
Windows; у Linux їх не повторювали.

## Стан

- Етап перекладу PR #89 завершено на `REVIEW_STOP_AAM_UA_TRANSLATION`.
- Редакторський прохід PR #90: див. [редакторський аудит](audit/PRE_CLAUDE_REVIEW.md).
- AAM-UA7: рисунки реконструйовано, типографіку гармонізовано, виконано незалежний мовний огляд. Стоп-точка
  `REVIEW_STOP_AAM_UA7`.

Авторами робочого рукопису за вказівкою дослідника зазначено Вячеслава Богданова та Данила Селіванова, у цьому
порядку. Наукові зміни та остаточне затвердження рукопису залишаються за дослідником.
