#include "gtest/gtest.h"

#include <memory>
#include <iostream>
#include <string>
#include <exception>

#include "../stream_evaluator.h"

class StreamEvaluatorTest : public ::testing::Test {
 protected:
  virtual void SetUp() {
    gold_iss = std::make_unique<std::istringstream>();
    temp_iss = std::make_unique<std::istringstream>();
  }

  std::unique_ptr<std::istringstream> gold_iss;
  std::unique_ptr<std::istringstream> temp_iss;

  btest::StreamEvaluator test_eval;
};

TEST_F(StreamEvaluatorTest, BadGoldStream)
{
  temp_iss->setstate(std::ios_base::goodbit);
  gold_iss->setstate(std::ios_base::badbit);
  test_eval.AdoptStreams(std::move(gold_iss), std::move(temp_iss));
  ASSERT_FALSE(test_eval.GoldGood());
  ASSERT_TRUE(test_eval.TempGood());
}

TEST_F(StreamEvaluatorTest, BadTempStream)
{
  gold_iss->setstate(std::ios_base::goodbit);
  temp_iss->setstate(std::ios_base::badbit);
  test_eval.AdoptStreams(std::move(gold_iss), std::move(temp_iss));
  ASSERT_FALSE(test_eval.TempGood());
  ASSERT_TRUE(test_eval.GoldGood());
}

TEST_F(StreamEvaluatorTest, BadBothStream)
{
  gold_iss->setstate(std::ios_base::badbit);
  temp_iss->setstate(std::ios_base::badbit);
  test_eval.AdoptStreams(std::move(gold_iss), std::move(temp_iss));
  ASSERT_FALSE(test_eval.TempGood());
  ASSERT_FALSE(test_eval.GoldGood());
}

TEST_F(StreamEvaluatorTest, SameStream)
{
  std::string input_text = "1\n2\n3\n4\n5";
  gold_iss->str(input_text);
  temp_iss->str(input_text);
  test_eval.AdoptStreams(std::move(gold_iss), std::move(temp_iss));
  ASSERT_TRUE(test_eval.Compare());
}

TEST_F(StreamEvaluatorTest, DiffStream)
{
  std::string gold_text = "1\n2\n3\n4\n5";
  std::string temp_text = "1\n2\nX\n4\n5";
  gold_iss->str(gold_text);
  temp_iss->str(temp_text);
  test_eval.AdoptStreams(std::move(gold_iss), std::move(temp_iss));
  ASSERT_FALSE(test_eval.Compare());
}

TEST_F(StreamEvaluatorTest, BadGoldStreamCompare)
{
  temp_iss->setstate(std::ios_base::goodbit);
  gold_iss->setstate(std::ios_base::badbit);
  test_eval.AdoptStreams(std::move(gold_iss), std::move(temp_iss));
  ASSERT_THROW(test_eval.Compare(), std::runtime_error);
}

TEST_F(StreamEvaluatorTest, BadTempStreamCompare)
{
  gold_iss->setstate(std::ios_base::goodbit);
  temp_iss->setstate(std::ios_base::badbit);
  test_eval.AdoptStreams(std::move(gold_iss), std::move(temp_iss));
  ASSERT_THROW(test_eval.Compare(), std::runtime_error);
}

TEST_F(StreamEvaluatorTest, CompareLongerGold)
{
  std::string gold_text = "1\n2\n3\n4\n5";
  std::string temp_text = "1\n2\n3\n4";
  gold_iss->str(gold_text);
  temp_iss->str(temp_text);
  test_eval.AdoptStreams(std::move(gold_iss), std::move(temp_iss));
  ASSERT_FALSE(test_eval.Compare());

}

TEST_F(StreamEvaluatorTest, CompareLongerTemp)
{
  std::string gold_text = "1\n2\n3\n4";
  std::string temp_text = "1\n2\n3\n4\n5";
  gold_iss->str(gold_text);
  temp_iss->str(temp_text);
  test_eval.AdoptStreams(std::move(gold_iss), std::move(temp_iss));
  ASSERT_FALSE(test_eval.Compare());
}

TEST_F(StreamEvaluatorTest, DiffWorks)
{
  std::string gold_text = "1\n2\nX\n4";
  std::string temp_text = "1\n2\n3\n4";
  std::string diff;
  gold_iss->str(gold_text);
  temp_iss->str(temp_text);

  test_eval.AdoptStreams(std::move(gold_iss), std::move(temp_iss));
  diff = test_eval.GetDiff();
  ASSERT_EQ(diff,"@@ -1,4 +1,4 @@\n 1\n 2\n-X\n+3\n 4\n");
}

TEST_F(StreamEvaluatorTest, BadGoldStreamDiff)
{
  temp_iss->setstate(std::ios_base::goodbit);
  gold_iss->setstate(std::ios_base::badbit);
  test_eval.AdoptStreams(std::move(gold_iss), std::move(temp_iss));
  ASSERT_THROW(test_eval.GetDiff(), std::runtime_error);
}

TEST_F(StreamEvaluatorTest, BadTempStreamDiff)
{
  gold_iss->setstate(std::ios_base::goodbit);
  temp_iss->setstate(std::ios_base::badbit);
  test_eval.AdoptStreams(std::move(gold_iss), std::move(temp_iss));
  ASSERT_THROW(test_eval.GetDiff(), std::runtime_error);
}
