/* Copyright (C) 2016 Philipp Benner
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

package gonetics

/* -------------------------------------------------------------------------- */

//import "fmt"
import "math"
import "testing"

/* -------------------------------------------------------------------------- */

func TestTF1(t *testing.T) {

  tf  := EmptyTFMatrix()
  err := tf.ImportMatrix("tf_test.table")

  if err != nil {
    t.Error("TestTF1 failed")
  }
  if math.Abs(tf.Get('a', 0) - -0.308122295362332) > 1e-8 {
    t.Error("TestTF1 failed")
  }
}

func TestTF2(t *testing.T) {

  tf  := EmptyTFMatrix()
  err := tf.ImportMatrix("tf_test.table")

  if err != nil {
    t.Error("TestTF2 failed")
  }

  seq := []byte("cacgtg")

  score1, _ := PWM{tf}.Score(seq, false)
  score2, _ := PWM{tf}.Score(seq, true)

  if score1 != score2 {
    t.Error("TestTF2 failed")
  }
}

func TestTF3(t *testing.T) {

  tf  := EmptyTFMatrix()
  err := tf.ImportMatrix("tf_test.table")

  if err != nil {
    t.Error("TestTF3 failed")
  }

  seq := []byte("cacgtgaaaccctttgg")

  score1 := PWM{tf}.MaxScore(seq, false)
  score2 := PWM{tf}.MaxScore(seq, true )

  if score1 != score2 {
    t.Error("TestTF3 failed")
  }
}
