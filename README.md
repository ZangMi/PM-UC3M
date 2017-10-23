<!-- README.md is generated from README.Rmd. Please edit that file -->
Notes for *Predictive Modeling*
===============================

[![Travis-CI Build Status](https://travis-ci.org/egarpor/SSS2-UC3M.svg?branch=master)](https://travis-ci.org/egarpor/PM-UC3M) [![License](https://img.shields.io/badge/license-CC_BY--NC--SA_4.0-blue.svg)](https://creativecommons.org/licenses/by-nc-sa/4.0/)

Course overview
---------------

These are the notes for *Predictive Modeling* for the course 2017/2018. The subject is part of the [MSc in Big Data Analytics](https://www.uc3m.es/ss/Satellite/Postgrado/en/Detalle/Estudio_C/1371210340413/1371219633369/Master_in_Big_Data_Analytics) from [Carlos III University of Madrid](http://www.uc3m.es/).

The notes are available at <https://bookdown.org/egarpor/PM-UC3M>.

Here is a broad view of the **syllabus**:

1.  [Introduction](https://bookdown.org/egarpor/PM-UC3M/intro.html)
2.  [Linear models I](https://bookdown.org/egarpor/PM-UC3M/lm-i.html)
3.  [Linear models II](https://bookdown.org/egarpor/PM-UC3M/lm-ii.html)
4.  [Generalized linear models](https://bookdown.org/egarpor/PM-UC3M/glm.html)
5.  [Generalized additive models](https://bookdown.org/egarpor/PM-UC3M/gam.html)
6.  [Nonparametric regression](https://bookdown.org/egarpor/PM-UC3M/npreg.html)
7.  [Regression trees](https://bookdown.org/egarpor/PM-UC3M/trees.html)

**Office hours**: Thursdays from 19:15 to 20:15, at the lab in which the lesson took place.

Course evaluation
-----------------

Evaluation is done by means of **four group reports**, to be done in groups of 3 (preferable) or 2 students, plus a presentation of one randomly-selected report. The final grade (in the scale 0-10) is:

    0.70 * reports_grade + 0.30 * individual_presentation_grade

The `reports_grade` is the weighted grade of the four reports about the following topics:

1.  Linear models I and II (weight: 40%; **deadline: 2017/12/14**).
2.  Generalized linear models (weight: 20%; **deadline: 2017/12/21**).
3.  Generalized additive models (weight: 20%; **deadline: 2018/01/04**).
4.  Nonparametric regression (weight: 20%; **deadline: 2018/01/25**).

### Group projects

#### Groups

It is up to you to form the groups based on your grade expectations, affinity, complementary skills, etc. Take into account that all the students in a group will be graded evenly for the `reports_grade`.

Communicate the group compositions by filling this [Google Sheet](https://docs.google.com/a/uc3m.es/spreadsheets/d/10zWtuhpAEtfZs7tuL9wEPdGsLVVYT-MKDLudvuhEFCA/edit?usp=sharing) (you need to log in with your UC3M account). An example: if students A, B, and C form the group number `4` in the third project, add a `4` to the the rows of A, B, and C in the `Project III` column. Keep the same numbers across columns if the composition of the group does not change. It is possible to switch group members from project to project, but be advised that this might have undesirable consequences: a person switching groups might have to present in two different randomly-chosen presentations.

#### Report structure

The report must analyze a real dataset of your choice using the statistical methodology that we have seen in the lessons and labs. The purpose is to demonstrate that you understand and know how to apply and *interpret* the studied statistical techniques in a real-case scenario. Simulation studies (*i.e.*, not real data but artificially generated) are also allowed if they are meant to illustrate interesting features or insisghts of the methods.

Compulsory format guidelines:

-   Structure: title, authors, abstract, first section, second section and so on. Like [this](http://cje.oxfordjournals.org/content/38/2/257.full.pdf+html). Do not use a cover.
-   Font size: 11 points.
-   Spacing: single space, single column.
-   Length: **10 pages at most**. This includes tables, graphics, code, equations, appendices, etc. Everything must fit in 10 pages. Want extra space? You can buy it at `-1.00` points per extra page!
-   Code: do not include *all* the code in the report. Instead of, just add the parts you want to highlight (or that you believe are particularly interesting) and describe what the code is doing.
-   Style: be concise, specific, and insightful. Make a clever use of the space.

TODO Mandatory report structure:

1.  **Abstract**. Provide a concise summary of the project. It must not exceed 250 words.
2.  **Introduction**. State what is the problem to be studied. Provide some context, the question(s) that you want to address, a motivation of its importance, references, etc. Must not be *very* extensive.
3.  **Methods**. Describe in detail the methods you are going to use, showing that you know what you are writing about. Provide background, useful interpretations, connections with other methods, alternative point of views, remarks, . Optionally, go beyond what was taught in the lessons, etc.
4.  **Statistical analysis**. Make use of some of the aforementioned statistical techniques, the ones that are more convenient to your particular case study. You can choose between covering several at a more superficial level, or one or two in more depth. Justify their adequacy and obtain analyses, explaining how you did it, in the form of plots and summaries. Provide a critical discussion about the outputs and give insights about them.
5.  **Conclusions**. Summary of what was addressed in the project and of the most important conclusions. Takeaway messages. The conclusions are not required to be spectacular, but *fair and honest* in terms of what you have discovered.
6.  **References**. Refer to the sources of information that you have employed (for the data, for information on the data, for the methodology, for the statistical analyses, etc).

#### Grading

The report grading will be performed according to the following breakdown:

-   **Originality** of the problem studied and data acquisition process (up to 1.5 points).
-   **Statistical analyses presented** and their depth (up to 3 points). At least two different techniques should be employed (simple and multiple linear regression count as different, but the use of other techniques as well is mostly encouraged). Graded depending on their adequacy to the problem studied and the evidence you demonstrate about your knowledge of them.
-   Accuracy of the **interpretation** of the analyses (up to 2 points). Graded depending on the detail and rigor of the insights you elaborate from the statistical analyses.
-   **Reproducibility** of the study (1.0 point). Awarded if the code for reproducing the study, as well as the data, is provided in a ready-to-use way (e.g. the outputs from `R Commander`'s report mode along with the data).
-   **Presentation** of the report (1.5 point). This involves the correct usage of English, the readability, the conciseness and the overall presentation quality.
-   **Excellence** (2 bonus points). Awarded for creativity of the analysis, use of advanced statistical tools, use of points briefly covered in lessons/labs, advanced insights into the methods, completeness of the report, use of advanced presentation tools, etc. Only awarded if the sum of regular points is above 7.5. TODO

#### Hand in

Send **one** email per group to <edgarcia@est-econ.uc3m.es> with the subject "\[PM - Report X - Group Y\]" (replace X and Y as convenience) and write the group members in the body of the email. Attach all relevant scripts with the report in a single `.zip` file. Ensure easy reproducibility of analysis (e.g. include the data you employed).

**Deadlines**: in general, two weeks after we cover in the lessons the report's topic. The deadlines finish at **23:59** and the list of days can be seen above. Recall that the last deadline, the 25th of January at 23:59, is quite close to the presentation day, the 30th of January.

Need extra time? Buy it at expenses of a penalization in your grade: `-0.04` points per extra hour! But only until the **25th of January**, reports received after that date will not be graded. The maximum grade according to the hours passed from the deadline is given in the graph below.

![](README/README-unnamed-chunk-2-1.png)

### Report presentation

The **30th of January, from 09:00 to 21:00, in room 0.B.08** each group must present one of the four projects, which will be randomly selected. On the 25th of January, a list with the presentation topics and presentation schedule will be public. Restrictions, if any, must be communicated and justified in advance. The presentation will determine the `individual_presentation_grade`.

The evaluation will be carried out as follows

-   Each group must present one report, with an strict **15-minutes limit** for the whole presentation. Tailor your presentation to that duration and break evenly the time among the group members. **All students must talk for at least 5 minutes.**
-   All group members must be present in the presentation. Failure to be at the presentation results in `individual_presentation_grade <- 0`.
-   The **presentation must follow tightly the report contents**. You just need to explain properly what you did in the group project without adding new material. Take this into account when preparing the reports.
-   In the presentation, you may want to highlight your main contributions and show your understanding of the methods.
-   **Questions will be asked during and after the presentation of each student.** Questions may be about the theory, implementation, or other matters related with the report. You are assumed to know and be able to defend all these topics. Together with the presentation, the accurateness, conciseness, and detail in the answers will determine the grade.

### Academic fraud

Evidences of academic fraud will have serious consequences, such as a zero grade for the whole group and the reporting of the fraud detection to the pertinent academic authorities. Academic fraud includes (but is not limited to) plagiarism, use of sources without proper credit, project outsourcing, and the use of external tutoring not mentioned explicitly.

Contributions
-------------

Contributions, reporting of typos, and feedback on the notes are very welcome. Either send an email to <edgarcia@est-econ.uc3m.es> or (preferably) fork the repository, make your changes and open a pull request. Give me a reason for writing your name in a list of contributors!

License
-------

All material in this repository is licensed under [CC BY-NC-SA 4.0](https://creativecommons.org/licenses/by-nc-sa/4.0/).
