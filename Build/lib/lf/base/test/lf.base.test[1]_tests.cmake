add_test( ForwardIteratorTest.nonConst /u/magina/Documents/lehrfempp/Build/lib/lf/base/test/lf.base.test [==[--gtest_filter=ForwardIteratorTest.nonConst]==] --gtest_also_run_disabled_tests)
set_tests_properties( ForwardIteratorTest.nonConst PROPERTIES WORKING_DIRECTORY /u/magina/Documents/lehrfempp/Build/lib/lf/base/test)
add_test( ForwardIterator.constEntities /u/magina/Documents/lehrfempp/Build/lib/lf/base/test/lf.base.test [==[--gtest_filter=ForwardIterator.constEntities]==] --gtest_also_run_disabled_tests)
set_tests_properties( ForwardIterator.constEntities PROPERTIES WORKING_DIRECTORY /u/magina/Documents/lehrfempp/Build/lib/lf/base/test)
add_test( ForwardIterator.DefaultConstructible /u/magina/Documents/lehrfempp/Build/lib/lf/base/test/lf.base.test [==[--gtest_filter=ForwardIterator.DefaultConstructible]==] --gtest_also_run_disabled_tests)
set_tests_properties( ForwardIterator.DefaultConstructible PROPERTIES WORKING_DIRECTORY /u/magina/Documents/lehrfempp/Build/lib/lf/base/test)
add_test( Range.useInForLoop /u/magina/Documents/lehrfempp/Build/lib/lf/base/test/lf.base.test [==[--gtest_filter=Range.useInForLoop]==] --gtest_also_run_disabled_tests)
set_tests_properties( Range.useInForLoop PROPERTIES WORKING_DIRECTORY /u/magina/Documents/lehrfempp/Build/lib/lf/base/test)
add_test( base_forwardRangeTest.initializerList /u/magina/Documents/lehrfempp/Build/lib/lf/base/test/lf.base.test [==[--gtest_filter=base_forwardRangeTest.initializerList]==] --gtest_also_run_disabled_tests)
set_tests_properties( base_forwardRangeTest.initializerList PROPERTIES WORKING_DIRECTORY /u/magina/Documents/lehrfempp/Build/lib/lf/base/test)
add_test( RandomAccessIterator.nonConst /u/magina/Documents/lehrfempp/Build/lib/lf/base/test/lf.base.test [==[--gtest_filter=RandomAccessIterator.nonConst]==] --gtest_also_run_disabled_tests)
set_tests_properties( RandomAccessIterator.nonConst PROPERTIES WORKING_DIRECTORY /u/magina/Documents/lehrfempp/Build/lib/lf/base/test)
add_test( RandomAccessIterator.InteractWithForwardIterators /u/magina/Documents/lehrfempp/Build/lib/lf/base/test/lf.base.test [==[--gtest_filter=RandomAccessIterator.InteractWithForwardIterators]==] --gtest_also_run_disabled_tests)
set_tests_properties( RandomAccessIterator.InteractWithForwardIterators PROPERTIES WORKING_DIRECTORY /u/magina/Documents/lehrfempp/Build/lib/lf/base/test)
add_test( RandomAccessIterator.RandomAccessFunctionality /u/magina/Documents/lehrfempp/Build/lib/lf/base/test/lf.base.test [==[--gtest_filter=RandomAccessIterator.RandomAccessFunctionality]==] --gtest_also_run_disabled_tests)
set_tests_properties( RandomAccessIterator.RandomAccessFunctionality PROPERTIES WORKING_DIRECTORY /u/magina/Documents/lehrfempp/Build/lib/lf/base/test)
add_test( RefEl.dimensionCorrect /u/magina/Documents/lehrfempp/Build/lib/lf/base/test/lf.base.test [==[--gtest_filter=RefEl.dimensionCorrect]==] --gtest_also_run_disabled_tests)
set_tests_properties( RefEl.dimensionCorrect PROPERTIES WORKING_DIRECTORY /u/magina/Documents/lehrfempp/Build/lib/lf/base/test)
set( lf.base.test_TESTS ForwardIteratorTest.nonConst ForwardIterator.constEntities ForwardIterator.DefaultConstructible Range.useInForLoop base_forwardRangeTest.initializerList RandomAccessIterator.nonConst RandomAccessIterator.InteractWithForwardIterators RandomAccessIterator.RandomAccessFunctionality RefEl.dimensionCorrect)
