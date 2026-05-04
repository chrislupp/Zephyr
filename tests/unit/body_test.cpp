/*
    Project Pegasus
    Version 0.1.0

    Copyright 2024-2026, Christopher A. Lupp

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published
    by the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/

#include "pegasus/geometry/body.h"

#include <gtest/gtest.h>

namespace pegasus {
namespace {

TEST(BodyTest, DefaultConstruction) {
    Body body;
    EXPECT_TRUE(body.GetName().empty());
    EXPECT_EQ(body.GetNumSections(), 0u);
    EXPECT_EQ(body.GetNodeCount(), 20);
    EXPECT_DOUBLE_EQ(body.GetNodeSpacing(), 1.0);
    EXPECT_FALSE(body.IsYDuplicated());
}

TEST(BodyTest, NamedConstruction) {
    Body body("fuselage");
    EXPECT_EQ(body.GetName(), "fuselage");
}

TEST(BodyTest, SetName) {
    Body body("old_name");
    body.SetName("new_name");
    EXPECT_EQ(body.GetName(), "new_name");
}

TEST(BodyTest, AddAndGetSections) {
    Body body("fuselage");

    BodySection s1;
    s1.center = Vec3(0.0, 0.0, 0.0);
    s1.width = 0.5;
    s1.height = 0.3;
    body.AddSection(s1);

    BodySection s2;
    s2.center = Vec3(1.0, 0.0, 0.0);
    s2.width = 0.8;
    s2.height = 0.6;
    body.AddSection(s2);

    EXPECT_EQ(body.GetNumSections(), 2u);

    const auto& sections = body.GetSections();
    EXPECT_EQ(sections.size(), 2u);
    EXPECT_DOUBLE_EQ(sections[0].width, 0.5);
    EXPECT_DOUBLE_EQ(sections[1].width, 0.8);
}

TEST(BodyTest, NodeCount) {
    Body body;
    body.SetNodeCount(30);
    EXPECT_EQ(body.GetNodeCount(), 30);
}

TEST(BodyTest, NodeSpacing) {
    Body body;
    body.SetNodeSpacing(2.5);
    EXPECT_DOUBLE_EQ(body.GetNodeSpacing(), 2.5);
}

TEST(BodyTest, Translation) {
    Body body;
    Vec3 t(1.0, 2.0, 3.0);
    body.SetTranslation(t);

    EXPECT_DOUBLE_EQ(body.GetTranslation().x(), 1.0);
    EXPECT_DOUBLE_EQ(body.GetTranslation().y(), 2.0);
    EXPECT_DOUBLE_EQ(body.GetTranslation().z(), 3.0);
}

TEST(BodyTest, Scale) {
    Body body;
    Vec3 s(2.0, 3.0, 4.0);
    body.SetScale(s);

    EXPECT_DOUBLE_EQ(body.GetScale().x(), 2.0);
    EXPECT_DOUBLE_EQ(body.GetScale().y(), 3.0);
    EXPECT_DOUBLE_EQ(body.GetScale().z(), 4.0);
}

TEST(BodyTest, DefaultScale) {
    Body body;
    EXPECT_DOUBLE_EQ(body.GetScale().x(), 1.0);
    EXPECT_DOUBLE_EQ(body.GetScale().y(), 1.0);
    EXPECT_DOUBLE_EQ(body.GetScale().z(), 1.0);
}

TEST(BodyTest, DefaultTranslation) {
    Body body;
    EXPECT_DOUBLE_EQ(body.GetTranslation().x(), 0.0);
    EXPECT_DOUBLE_EQ(body.GetTranslation().y(), 0.0);
    EXPECT_DOUBLE_EQ(body.GetTranslation().z(), 0.0);
}

TEST(BodyTest, YDuplicate) {
    Body body;
    EXPECT_FALSE(body.IsYDuplicated());
    body.SetYDuplicate(true);
    EXPECT_TRUE(body.IsYDuplicated());
    body.SetYDuplicate(false);
    EXPECT_FALSE(body.IsYDuplicated());
}

}  // namespace
}  // namespace pegasus
