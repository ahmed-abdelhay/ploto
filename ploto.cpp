#include <assert.h>
#include <math.h>
#include <stdint.h>

#include <charconv>
#include <string>
#include <vector>

#include <GL/gl3w.h>
#include <GLFW/glfw3.h>
#include <imgui.h>
#include <imgui_impl_glfw.h>
#include <imgui_impl_opengl3.h>
#include <implot.h>

#define RETURN_IF_FALSE(expr)                                                  \
  do {                                                                         \
    if (!(expr)) {                                                             \
      return false;                                                            \
    }                                                                          \
  } while (0)

static const size_t MAX_EXPRESSION_SIZE = 512;
static const ImVec4 BLUE(41 / 255., 74 / 255., 122 / 255., 1);
static const ImVec4 RED(1, 0, 0, 1);

//---------------------------Parsing-------------------------
/*
E --> E + E
    | E - E
    | - E
    | E * E
    | E ** E
    | E / E
    | sin( E )
    | cos( E )
    | tan( E )
    | pow( E , E )
    | var
    | const
*/

struct Buffer {
  const char *data = nullptr;
  size_t size = 0;
  size_t cursor = 0;

  char Peek() const {
    assert(cursor < size);
    return data[cursor];
  }
};

enum class ExpTokenType {
  IDENTIFIER,
  NUMERIC_LITERAL,
  // operators
  OPERATOR_PLUS,
  OPERATOR_MINUS,
  OPERATOR_MULTIPLY,
  OPERATOR_DIVIDE,
  OPERATOR_POW,
  // functions
  FUNCTION_SIN,
  FUNCTION_COS,
  FUNCTION_TAN,
  FUNCTION_POW,
  //
  LEFT_PARAN,
  RIGHT_PARAN,
  COMMA,
};

struct Token {
  ExpTokenType type;
  // data.
  std::string name;      // when type == IDENTIFIER.
  double numericLiteral; // when type == NUMERIC_LITERAL.
};

struct Expression {
  enum class NodeType {
    ADD,
    SUB,
    MULTIPLY,
    POWER,
    DIVIDE,
    SIN,
    COS,
    TAN,
    VARIABLE,
    CONSTANAT
  };

  struct Node {
    NodeType type;
    int32_t children[2] = {-1, -1};
    std::string variable;
    double constant;
  };

  std::vector<Node> nodes;
  size_t rootNodeId = 0;

  double Evaluate(int32_t nodeId, double x) const {
    const Node &node = nodes[nodeId];
    switch (node.type) {
    case NodeType::ADD: {
      const double v0 = Evaluate(node.children[0], x);
      const double v1 = Evaluate(node.children[1], x);
      return v0 + v1;
    } break;
    case NodeType::SUB: {
      const double v0 = Evaluate(node.children[0], x);
      const double v1 = Evaluate(node.children[1], x);
      return v0 - v1;
    } break;
    case NodeType::DIVIDE: {
      const double v0 = Evaluate(node.children[0], x);
      const double v1 = Evaluate(node.children[1], x);
      return v0 / v1;
    } break;
    case NodeType::MULTIPLY: {
      const double v0 = Evaluate(node.children[0], x);
      const double v1 = Evaluate(node.children[1], x);
      return v0 * v1;
    } break;
    case NodeType::POWER: {
      const double v0 = Evaluate(node.children[0], x);
      const double v1 = Evaluate(node.children[1], x);
      return pow(v0, v1);
    } break;
    case NodeType::SIN: {
      const double v0 = Evaluate(node.children[0], x);
      return sin(v0);
    } break;
    case NodeType::COS: {
      const double v0 = Evaluate(node.children[0], x);
      return cos(v0);
    } break;
    case NodeType::TAN: {
      const double v0 = Evaluate(node.children[0], x);
      return tan(v0);
    } break;
    case NodeType::VARIABLE: {
      return x;
    } break;
    case NodeType::CONSTANAT: {
      return node.constant;
    } break;
    }
    assert(false);
    return 0;
  }
  double Evaluate(double x) const {
    assert(!nodes.empty());
    return Evaluate(rootNodeId, x);
  }
};

// lexer
enum class LexerErrorType { OK, PARSE_ERROR, MULTIPLE_VARIABLES };

struct Tokens {
  std::vector<Token> tokens;
  size_t cursor = 0;

  bool Next(Token &value) const {
    if (cursor < tokens.size()) {
      value = tokens[cursor];
      return true;
    }
    return false;
  }

  bool Consume() {
    if (cursor < tokens.size()) {
      cursor++;
      return true;
    }
    return false;
  }

  bool Expect(ExpTokenType type) {
    Token next;
    if (Next(next) && next.type == type) {
      Consume();
      return true;
    }
    return false;
  }
};

struct LexerResult {
  Tokens tokens;
  LexerErrorType error = LexerErrorType::PARSE_ERROR;
  size_t errorLocation = 0;
};

static bool CompareWordAndSkip(Buffer &buffer, const std::string_view word) {
  const size_t length = word.size();
  if (buffer.cursor + length <= buffer.size) {
    for (size_t i = 0; i < length; ++i) {
      if (buffer.data[buffer.cursor + i] != word[i]) {
        return false;
      }
    }
    buffer.cursor += length;
    return true;
  }
  return false;
}

static bool ParseFloat(Buffer &buffer, double &result) {
  if (isdigit(buffer.Peek())) {
    std::string data;
    bool dotAdded = false;
    while (buffer.cursor < buffer.size &&
           (isdigit(buffer.Peek()) || buffer.Peek() == '.')) {
      if (buffer.Peek() == '.') {
        if (dotAdded) {
          return false;
        }
        dotAdded = true;
      }
      data.push_back(buffer.Peek());
      buffer.cursor++;
    }
    return std::from_chars(data.data(), data.data() + data.size(), result).ec ==
           std::errc();
  }
  return false;
}

static bool ParseIdentifier(Buffer &buffer, std::string &result) {
  if (isalpha(buffer.Peek())) {
    result.clear();
    while (buffer.cursor < buffer.size &&
           (isalnum(buffer.Peek()) || buffer.Peek() == '_')) {
      result.push_back(buffer.Peek());
      buffer.cursor++;
    }
    return true;
  }
  return false;
}

static LexerResult Tokenize(Buffer &buffer) {
  struct TokenString {
    ExpTokenType type;
    const char *string = nullptr;
  };

  static const TokenString tokensStrings[]{
      {ExpTokenType::FUNCTION_SIN, "sin"},
      {ExpTokenType::FUNCTION_COS, "cos"},
      {ExpTokenType::FUNCTION_TAN, "tan"},
      {ExpTokenType::FUNCTION_POW, "pow"},
      {ExpTokenType::OPERATOR_PLUS, "+"},
      {ExpTokenType::OPERATOR_MINUS, "-"},
      {ExpTokenType::OPERATOR_POW, "**"},
      {ExpTokenType::OPERATOR_MULTIPLY, "*"},
      {ExpTokenType::OPERATOR_DIVIDE, "/"},
      {ExpTokenType::LEFT_PARAN, "("},
      {ExpTokenType::RIGHT_PARAN, ")"},
      {ExpTokenType::COMMA, ","},
  };

  constexpr size_t count = std::size(tokensStrings);
  LexerResult result;
  while (buffer.cursor < buffer.size) {
    Token token;
    while (buffer.cursor < buffer.size && isblank(buffer.Peek())) {
      buffer.cursor++;
    }
    if (buffer.cursor >= buffer.size) {
      continue;
    }

    bool found = false;
    for (size_t i = 0; i < count; ++i) {
      if (CompareWordAndSkip(buffer, tokensStrings[i].string)) {
        token.type = tokensStrings[i].type;
        result.tokens.tokens.push_back(token);
        found = true;
        break;
      }
    }
    if (found) {
      continue;
    }
    if (ParseIdentifier(buffer, token.name)) {
      token.type = ExpTokenType::IDENTIFIER;
      result.tokens.tokens.push_back(token);
    } else if (ParseFloat(buffer, token.numericLiteral)) {
      token.type = ExpTokenType::NUMERIC_LITERAL;
      result.tokens.tokens.push_back(token);
    } else {
      result.error = LexerErrorType::PARSE_ERROR;
      result.errorLocation = buffer.cursor;
      return result;
    }
  }

  // specific to this application: only accept at most a single variable
  // functions and the variable name is x
  for (const Token token : result.tokens.tokens) {
    if (token.type == ExpTokenType::IDENTIFIER && token.name != "x") {
      result.error = LexerErrorType::MULTIPLE_VARIABLES;
      return result;
    }
  }

  result.error = LexerErrorType::OK;
  return result;
}

// parser
static bool IsBinaryOperator(ExpTokenType type) {
  return type == ExpTokenType::OPERATOR_MULTIPLY ||
         type == ExpTokenType::OPERATOR_DIVIDE ||
         type == ExpTokenType::OPERATOR_PLUS ||
         type == ExpTokenType::OPERATOR_POW ||
         type == ExpTokenType::OPERATOR_MINUS;
}
static bool ParseP(Tokens &tokens, Expression &result);
static bool ParseE(Tokens &tokens, Expression &result);

static bool ParseE(Tokens &tokens, Expression &result) {
  Expression::Node node;
  node.children[0] = result.nodes.size();
  RETURN_IF_FALSE(ParseP(tokens, result));
  Token next;
  while (tokens.Next(next) && IsBinaryOperator(next.type)) {
    RETURN_IF_FALSE(tokens.Consume());
    if (next.type == ExpTokenType::OPERATOR_MULTIPLY) {
      node.type = Expression::NodeType::MULTIPLY;
    } else if (next.type == ExpTokenType::OPERATOR_PLUS) {
      node.type = Expression::NodeType::ADD;
    } else if (next.type == ExpTokenType::OPERATOR_DIVIDE) {
      node.type = Expression::NodeType::DIVIDE;
    } else if (next.type == ExpTokenType::OPERATOR_MINUS) {
      node.type = Expression::NodeType::SUB;
    } else if (next.type == ExpTokenType::OPERATOR_POW) {
      node.type = Expression::NodeType::POWER;
    }
    node.children[1] = result.nodes.size();
    RETURN_IF_FALSE(ParseP(tokens, result));
    result.rootNodeId = result.nodes.size();
    result.nodes.push_back(node);
  }
  return true;
}

static bool ParseP(Tokens &tokens, Expression &result) {
  Token next;
  if (tokens.Next(next) && next.type == ExpTokenType::IDENTIFIER) {
    RETURN_IF_FALSE(tokens.Consume());
    Expression::Node node;
    node.variable = next.name;
    node.type = Expression::NodeType::VARIABLE;
    result.nodes.push_back(node);
    return true;
  } else if (tokens.Next(next) && next.type == ExpTokenType::NUMERIC_LITERAL) {
    RETURN_IF_FALSE(tokens.Consume());
    Expression::Node node;
    node.constant = next.numericLiteral;
    node.type = Expression::NodeType::CONSTANAT;
    result.nodes.push_back(node);
    return true;
  } else if (tokens.Next(next) && next.type == ExpTokenType::FUNCTION_COS) {
    RETURN_IF_FALSE(tokens.Consume());
    RETURN_IF_FALSE(tokens.Expect(ExpTokenType::LEFT_PARAN));
    const size_t nodeId = result.nodes.size();
    result.nodes.push_back({});
    result.nodes[nodeId].type = Expression::NodeType::COS;
    result.nodes[nodeId].children[0] = result.nodes.size();
    RETURN_IF_FALSE(ParseE(tokens, result));
    RETURN_IF_FALSE(tokens.Expect(ExpTokenType::RIGHT_PARAN));
    return true;
  } else if (tokens.Next(next) && next.type == ExpTokenType::FUNCTION_SIN) {
    RETURN_IF_FALSE(tokens.Consume());
    RETURN_IF_FALSE(tokens.Expect(ExpTokenType::LEFT_PARAN));
    const size_t nodeId = result.nodes.size();
    result.nodes.push_back({});
    result.nodes[nodeId].type = Expression::NodeType::SIN;
    result.nodes[nodeId].children[0] = result.nodes.size();
    RETURN_IF_FALSE(ParseE(tokens, result));
    RETURN_IF_FALSE(tokens.Expect(ExpTokenType::RIGHT_PARAN));
    return true;
  } else if (tokens.Next(next) && next.type == ExpTokenType::FUNCTION_TAN) {
    RETURN_IF_FALSE(tokens.Consume());
    RETURN_IF_FALSE(tokens.Expect(ExpTokenType::LEFT_PARAN));
    const size_t nodeId = result.nodes.size();
    result.nodes.push_back({});
    result.nodes[nodeId].type = Expression::NodeType::TAN;
    result.nodes[nodeId].children[0] = result.nodes.size();
    RETURN_IF_FALSE(ParseE(tokens, result));
    RETURN_IF_FALSE(tokens.Expect(ExpTokenType::RIGHT_PARAN));
    return true;
  } else if (tokens.Next(next) && next.type == ExpTokenType::FUNCTION_POW) {
    RETURN_IF_FALSE(tokens.Consume());
    RETURN_IF_FALSE(tokens.Expect(ExpTokenType::LEFT_PARAN));
    const size_t nodeId = result.nodes.size();
    result.nodes.push_back({});
    result.nodes[nodeId].type = Expression::NodeType::POWER;
    result.nodes[nodeId].children[0] = result.nodes.size();
    RETURN_IF_FALSE(ParseE(tokens, result));
    RETURN_IF_FALSE(tokens.Expect(ExpTokenType::COMMA));
    result.nodes[nodeId].children[1] = result.nodes.size();
    RETURN_IF_FALSE(ParseE(tokens, result));
    RETURN_IF_FALSE(tokens.Expect(ExpTokenType::RIGHT_PARAN));
    return true;
  } else if (tokens.Next(next) && next.type == ExpTokenType::LEFT_PARAN) {
    RETURN_IF_FALSE(tokens.Consume());
    RETURN_IF_FALSE(ParseE(tokens, result));
    RETURN_IF_FALSE(tokens.Expect(ExpTokenType::RIGHT_PARAN));
    return true;
  }
  return false;
}

static bool ParseTokens(Tokens tokens, Expression &result) {
  if (ParseE(tokens, result)) {
    Token next;
    return !tokens.Next(next);
  }
  return false;
}

/*
 * TODO:
 *   - Write the actual parsing code.
 *   - Add controls for the type of plot and the range of the plot in x and y
 * directions.
 */
struct Color {
  uint8_t r = 0;
  uint8_t g = 0;
  uint8_t b = 0;
  uint8_t a = 255;
};

static Color GenerateColor() {
  static Color pregenerateList[]{
      Color{255, 0, 0, 255},   Color{255, 255, 0, 255}, Color{255, 0, 255, 255},
      Color{0, 255, 255, 255}, Color{0, 0, 255, 255},   Color{0, 255, 0, 255},
  };

  static const size_t listSize = sizeof(pregenerateList) / sizeof(Color);
  static size_t counter = 0;

  const Color result = pregenerateList[counter];
  counter = (counter + 1) % listSize;
  return result;
}

struct Function {
  Expression expr;
  char name[MAX_EXPRESSION_SIZE] = {};
  Color color = GenerateColor();
  bool visible = true;
  int bounds[2] = {0, 1000};
  float sampleDistance = 0.1;
};

static bool GenerateExpression(const char *text, Expression &result,
                               std::string &errorLog) {
  Buffer buffer;
  buffer.data = text;
  buffer.size = strlen(text);
  const LexerResult lexerResults = Tokenize(buffer);
  switch (lexerResults.error) {
  case LexerErrorType::MULTIPLE_VARIABLES:
    errorLog = "Only one functions with 1 variables are supported";
    return false;
    break;
  case LexerErrorType::PARSE_ERROR:
    errorLog = "Error parsing the expression at index:" +
               std::to_string(lexerResults.errorLocation);
    return false;
    break;
  case LexerErrorType::OK:
    break;
  }
  if (!ParseTokens(lexerResults.tokens, result)) {
    errorLog = "Error creating the expression tree";
    return false;
  }
  return true;
}

static void SampleFunction(const Function &f, std::vector<float> &x,
                           std::vector<float> &y) {
  const size_t samples = ((f.bounds[1] - f.bounds[0]) + 1) / f.sampleDistance;
  x.clear();
  y.clear();
  x.reserve(samples);
  y.reserve(samples);
  for (float t = f.bounds[0]; t <= f.bounds[1]; t += f.sampleDistance) {
    x.push_back(t);
    y.push_back(f.expr.Evaluate(t));
  }
}

struct State {
  struct GuiState {
    bool error = false;
    std::string errorLog;
    Function function;
  };

  struct PlotState {
    std::vector<float> x, y;
  };

  std::vector<Function> functions;
  std::vector<PlotState> plots;
  GuiState gui;
  bool redraw = true;

  void Update() {
    ImGuiStyle &style = ImGui::GetStyle();
    style.FrameRounding = style.GrabRounding = 12;

    ImGui::SetNextWindowPos(ImVec2(0, 0));
    ImGui::SetNextWindowSize(ImGui::GetIO().DisplaySize);

    ImGui::Begin("Viewer", nullptr,
                 ImGuiWindowFlags_NoTitleBar | ImGuiWindowFlags_NoScrollbar |
                     ImGuiWindowFlags_NoMove | ImGuiWindowFlags_MenuBar |
                     ImGuiWindowFlags_NoBringToFrontOnFocus |
                     ImGuiWindowFlags_NoResize);
    const ImVec2 minPoint = ImGui::GetWindowContentRegionMin();
    const ImVec2 maxPoint = ImGui::GetWindowContentRegionMax();
    const float height = (maxPoint.y - minPoint.y);

    ImGui::BeginChild("Plot view", ImVec2(-1, height * 0.8));
    if (ImPlot::BeginPlot("Plot", ImVec2(-1, height * 0.8))) {
      if (redraw) {
        plots.clear();
        for (const Function &f : functions) {
          PlotState p;
          SampleFunction(f, p.x, p.y);
          plots.push_back(p);
        }
        redraw = false;
      }
      assert(functions.size() == plots.size());
      for (size_t i = 0; i < plots.size(); i++) {
        if (functions[i].visible) {
          const ImVec4 color(functions[i].color.r / 255.0f,
                             functions[i].color.g / 255.0f,
                             functions[i].color.b / 255.0f, 1);
          ImPlot::PushStyleColor(ImPlotCol_Line, color);
          ImPlot::PlotLine(functions[i].name, plots[i].x.data(),
                           plots[i].y.data(), plots[i].x.size());
          ImPlot::PopStyleColor();
        }
      }
      ImPlot::EndPlot();
    }
    ImGui::EndChild();
    ImGui::BeginChild("Controls", ImVec2(-1, height * 0.2));
    {
      ImGui::Text("Function:");
      ImGui::SameLine();
      ImGui::InputText("##", gui.function.name, MAX_EXPRESSION_SIZE);
      ImGui::Text("X min:");
      ImGui::SameLine();
      ImGui::DragInt("##0", &gui.function.bounds[0], 1, -1000, 1000);
      ImGui::Text("X max:");
      ImGui::SameLine();
      ImGui::DragInt("##1", &gui.function.bounds[1], 1, -1000, 1000);
      ImGui::Text("Sampling distance:");
      ImGui::SameLine();
      ImGui::DragFloat("##2", &gui.function.sampleDistance, 1, 0, 10);
      if (ImGui::Button("Plot")) {
        if (GenerateExpression(gui.function.name, gui.function.expr,
                               gui.errorLog)) {
          functions.push_back(gui.function);
          redraw = true;
          gui = {};
        } else {
          gui.error = true;
        }
      }

      if (gui.error) {
        ImGui::TextColored(RED, gui.errorLog.c_str());
      }

      ImGui::TextColored(BLUE, "Functions");
      for (Function &function : functions) {
        ImGui::Spacing();
        if (ImGui::Checkbox(function.name, &function.visible)) {
          redraw = true;
        }
        ImGui::SameLine();
        ImGui::PushID(function.name);
        float color[3] = {function.color.r / 255.0f, function.color.g / 255.0f,
                          function.color.b / 255.0f};
        if (ImGui::ColorEdit3("Colour", color,
                              ImGuiColorEditFlags_NoInputs |
                                  ImGuiColorEditFlags_NoLabel)) {
          function.color.r = color[0] * 255;
          function.color.g = color[1] * 255;
          function.color.b = color[2] * 255;
          redraw = true;
        }
        ImGui::PopID();
      }
    }
    ImGui::Text("Application average: %.1f FPS", ImGui::GetIO().Framerate);
    ImGui::EndChild();
    ImGui::End();
  }
};

int main(int argc, char **argv) {
  // Setup window
  if (!glfwInit()) {
    return EXIT_FAILURE;
  }
  // GL 3.0 + GLSL 130
  const char *glsl_version = "#version 130";
  glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
  glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 0);

  // Create window with graphics context
  GLFWwindow *window = glfwCreateWindow(1280, 720, "Ploto", NULL, NULL);
  if (!window) {
    return EXIT_FAILURE;
  }
  glfwMakeContextCurrent(window);
  glfwSwapInterval(1); // Enable vsync

  // Initialize OpenGL loader
  if (gl3wInit() != 0) {
    return EXIT_FAILURE;
  }

  // Setup Dear ImGui context
  IMGUI_CHECKVERSION();
  ImGui::CreateContext();
  ImPlot::CreateContext();
  ImGuiIO &io = ImGui::GetIO();

  // Setup Dear ImGui style
  ImGui::StyleColorsDark();

  // Setup Platform/Renderer backends
  ImGui_ImplGlfw_InitForOpenGL(window, true);
  ImGui_ImplOpenGL3_Init(glsl_version);

  // Load Fonts
  if (ImFont *font = io.Fonts->AddFontFromFileTTF("d:/Fira.ttf", 25.0f)) {
    io.FontDefault = font;
  }

  // Our state
  static const ImVec4 clearColor(0.45f, 0.55f, 0.60f, 1.00f);
  State state;
  // Main loop
  while (!glfwWindowShouldClose(window)) {
    glfwPollEvents();

    // Start the Dear ImGui frame
    ImGui_ImplOpenGL3_NewFrame();
    ImGui_ImplGlfw_NewFrame();
    ImGui::NewFrame();

    state.Update();

    // Rendering
    ImGui::Render();
    int width, height;
    glfwGetFramebufferSize(window, &width, &height);
    glViewport(0, 0, width, height);
    glClearColor(clearColor.x, clearColor.y, clearColor.z, clearColor.w);
    glClear(GL_COLOR_BUFFER_BIT);
    ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
    glfwSwapBuffers(window);
  }

  // Cleanup
  ImGui_ImplOpenGL3_Shutdown();
  ImGui_ImplGlfw_Shutdown();
  ImPlot::DestroyContext();
  ImGui::DestroyContext();

  glfwDestroyWindow(window);
  glfwTerminate();

  return EXIT_SUCCESS;
}
