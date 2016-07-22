import hddm_x
import sys

def report(xrec):
   courses = xrec.getCourses()
   total_courses = len(courses)
   total_enrolled = 0
   total_credits = 0
   total_passed = 0
   for course in courses:
      if course.getResult().Pass:
         if course.year > 1992:
            total_credits += course.credits
         total_passed += 1
   print courses[0].name, "enrolled in", total_courses, "courses", \
         "and passed", total_passed, "of them,"
   print "earning a total of", total_credits, "credits."

ostr = hddm_x.ostream("exam2.hddm")
#ostr.compression = hddm_x.k_z_compression

for n in range(0, 1000000):
   xrec = hddm_x.HDDM()
   student = xrec.addStudents()
   student[0].name = "Humphrey Gaston"
   enrolled = student[0].addEnrolleds()
   enrolled[0].year = 2005
   enrolled[0].semester = 2
   course = enrolled[0].addCourses(3)
   course[0].credits = 3
   course[0].title = "Beginning Russian"
   course[0].addResults()
   course[0].getResult().grade = "A-"
   course[0].getResult().Pass = True
   course[1].credits = 1
   course[1].title = "Bohemian Poetry"
   course[1].addResults()
   course[1].getResult().grade = "C"
   course[1].getResult().Pass = 1
   course[2].credits = 4
   course[2].title = "Developmental Psychology"
   course[2].addResults()
   course[2].getResult().grade = "B+"
   course[2].getResult().Pass = True
   ostr.write(xrec)
del ostr

# try reading it back in from the output file just written

istr = hddm_x.istream("exam2.hddm")
count=0
for rec in istr:
   if count % 100000 == 0:
      print "event", count
      report(rec)
   count += 1
print "finished after", count, "events read."
