#include <hddm_x.h>

int report(x_HDDM_t *tscript);

int main()
{
   x_iostream_t* fp;
   x_HDDM_t* tscript;
   x_Student_t* student;
   x_Enrolleds_t* enrolleds;
   x_Courses_t* courses;
   x_Result_t* result;
   string_t name;
   string_t grade;
   string_t course;
   int i;
   
   fp = init_x_HDDM("exam2.hddm");
   for (i=0;i<1000000;i++) {
   
      // first build the complete nodal structure for this record
      tscript = make_x_HDDM();
      tscript->student = student = make_x_Student();
      student->enrolleds = enrolleds = make_x_Enrolleds(99);
      enrolleds->mult = 1;
      enrolleds->in[0].courses = courses = make_x_Courses(99);
      courses->mult = 3;
      courses->in[0].result = make_x_Result();
      courses->in[1].result = make_x_Result();
      courses->in[2].result = make_x_Result();
      
      // now fill in the attribute data for this record
      name = (char*)malloc(30);
      strcpy(name,"Humphrey Gaston");
      student->name = name;
      enrolleds->in[0].year = 2005;
      enrolleds->in[0].semester = 2;
      courses->in[0].credits = 3;
      course = (char*)malloc(30);
      courses->in[0].title = strcpy(course,"Beginning Russian");
      grade = (char*)malloc(5);
      courses->in[0].result->grade = strcpy(grade,"A-");
      courses->in[0].result->Pass = 1;
      courses->in[1].credits = 1;
      course = (char*)malloc(30);
      courses->in[1].title = strcpy(course,"Bohemian Poetry");
      grade = (char*)malloc(5);
      courses->in[1].result->grade = strcpy(grade,"C");
      courses->in[1].result->Pass = 1;
      courses->in[2].credits = 4;
      course = (char*)malloc(30);
      courses->in[2].title = strcpy(course,"Developmental Psychology");
      grade = (char*)malloc(5);
      courses->in[2].result->grade = strcpy(grade,"B+");
      courses->in[2].result->Pass = 1;
      
      // now open a file and write this one record into it
      flush_x_HDDM(tscript,fp);
   }
   close_x_HDDM(fp);

   // try reading it back in from the output file just written
   fp = open_x_HDDM("exam2.hddm");
   int count=0;
   while ((tscript = read_x_HDDM(fp)) != 0) {
      if (count/100000*100000 == count) {
         printf("event %d\n", count);
         report(tscript);
      }
      ++count;
   }
   printf("finished after %d events read.\n", count);
   return 0;
}

int report(x_HDDM_t *tscript)
{
   x_Student_t* student = tscript->student;
   x_Enrolleds_t* enrolleds = student->enrolleds;
   x_Courses_t* courses = enrolleds->in[0].courses;
   int total_courses = courses->mult;
   int total_enrolled = 0;
   int total_credits = 0;
   int total_passed = 0;
   int i;
   for (i=0;i<total_courses;++i) {
      if (courses->in[i].result->Pass) {
         if (enrolleds->in[0].year > 1992) {
            total_credits += courses->in[i].credits;
         }
         ++total_passed;
      }
   }
   printf("%s enrolled in %d courses and passed %d of them,\n",
          student->name, total_courses, total_passed);
   printf("earning a total of %d credits.\n", total_credits);
   return 0;
}
