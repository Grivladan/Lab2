#pragma once

namespace flow_simulator {

	using namespace System;
	using namespace System::ComponentModel;
	using namespace System::Collections;
	using namespace System::Windows::Forms;
	using namespace System::Data;
	using namespace System::Drawing;
	using namespace System::Drawing::Drawing2D;

	/// <summary>
	/// Summary for MainForm
	///
	/// WARNING: If you change the name of this class, you will need to change the
	///          'Resource File Name' property for the managed resource compiler tool
	///          associated with all .resx files this class depends on.  Otherwise,
	///          the designers will not be able to interact properly with localized
	///          resources associated with this form.
	/// </summary>
	public ref class MainForm : public System::Windows::Forms::Form
	{
	private:
		BufferedGraphics^ bufferedGraphics;
		BufferedGraphicsContext^ bufferedGraphicsContext;
		GraphicsState^ graphicsState;

		float alpha;
		int M;
		float gamma0;

		float zoom;
		Vector2f* offset;
		bool isDragged;
		int prevX, prevY;
		int step;

		float minValue, maxValue;
		std::vector< std::vector< float > >* value;

		FlowSolver* solver;

	private: System::Windows::Forms::ToolTip^  toolTip;
	private: System::Windows::Forms::NumericUpDown^  alphaInputBox;
	private: System::Windows::Forms::Label^  label1;
	private: System::Windows::Forms::NumericUpDown^  MInputBox;
	private: System::Windows::Forms::NumericUpDown^  Gamma0InputBox;
	private: System::Windows::Forms::Label^  label2;
	private: System::Windows::Forms::Label^  label3;

	private: System::Windows::Forms::Button^  drawPhiButton;
	private: System::Windows::Forms::Button^  drawPsiButton;
	private: System::Windows::Forms::Button^  drawCpButton;
	private: System::Windows::Forms::Button^  flowOnlyButton;

	private: System::Windows::Forms::Label^  label5;
	private: System::Windows::Forms::Label^  label6;
	private: System::Windows::Forms::Label^  label7;
	private: System::Windows::Forms::Label^  label8;
	private: System::Windows::Forms::Label^  label9;
	private: System::Windows::Forms::Label^  label10;
	private: System::Windows::Forms::Label^  phiLabel;
	private: System::Windows::Forms::Label^  psiLabel;
	private: System::Windows::Forms::Label^  VLabel;
	private: System::Windows::Forms::Label^  CpLabel;
	private: System::Windows::Forms::Label^  xLabel;
	private: System::Windows::Forms::Label^  yLabel;
	private: System::Windows::Forms::GroupBox^  ControlsGroupBox;

	private: System::Windows::Forms::Label^  VScalLabel;
	private: System::Windows::Forms::Label^  label12;
	private: System::Windows::Forms::Button^  VScalButton;
	private: System::Windows::Forms::Label^  YMaxLabel;
	private: System::Windows::Forms::Label^  XMinLabel;
	private: cli::array<Windows::Forms::Label^, 1>^ GradientLabels;

	private: System::Windows::Forms::Label^  YMinLabel;

	private: System::Windows::Forms::Label^  XMaxLabel;
	private: System::Windows::Forms::Label^  YMidLabel;
	private: System::Windows::Forms::Label^  XMidLabel;

	private: System::Windows::Forms::Timer^  animationTimer;
	private: System::Windows::Forms::Label^  timeLabel;
	private: System::Windows::Forms::Label^  label11;
	private: System::Windows::Forms::Label^  label13;
	private: System::Windows::Forms::NumericUpDown^  stepsBox;
	private: System::Windows::Forms::Button^  AutoStepsButton;
	private: System::Windows::Forms::Button^  RestartButton;
	private: System::Windows::Forms::Button^  NextStepButton;

	private: System::Windows::Forms::Timer^  starterTimer;

	public:
		MainForm(void)
		{
			InitializeComponent();

			bufferedGraphicsContext = BufferedGraphicsManager::Current;
			bufferedGraphicsContext->MaximumBuffer = Drawing::Size(ImageBox->Width, ImageBox->Height);
			bufferedGraphics = bufferedGraphicsContext->Allocate(ImageBox->CreateGraphics(), 
				System::Drawing::Rectangle(Point(0, 0), Drawing::Size(ImageBox->Width, ImageBox->Height)));

			alpha = 0.0f;
			M = 61;
			gamma0 = 1.0f;

			zoom = 300.0f;
			offset = new Vector2f(0.0f, 0.0f);
			*offset = pixelToReal(std::make_pair(-ImageBox->Width / 2, 3 * ImageBox->Height / 2));
			isDragged = false;

			solver = new FlowSolver();
			solver->setAlpha(alpha);
			solver->setGamma0(gamma0);
			solver->setM(M);
			solver->init();

			minValue = Inf;
			maxValue = - Inf;
			value = new std::vector< std::vector< float > >;

			updateAxisLabels();

			GradientLabels = gcnew cli::array<Windows::Forms::Label^, 1>(25);
			for (int i = 0; i < 25; i++)
			{
				Label^ t = gcnew System::Windows::Forms::Label();
				t->AutoSize = true;
				t->Location = System::Drawing::Point(
					ImageBox->Right + 2,
					ImageBox->Top - 6 + i * 560 / 24);
				t->TabIndex = 5;
				t->Anchor = static_cast<System::Windows::Forms::AnchorStyles>(System::Windows::Forms::AnchorStyles::Top 
					| System::Windows::Forms::AnchorStyles::Right);
				this->ImageGroupBox->Controls->Add(t);
				GradientLabels[i] = t;
			}

			redraw(0);
		}

	protected:
		/// <summary>
		/// Clean up any resources being used.
		/// </summary>
		~MainForm()
		{
			if (components)
			{
				delete components;
			}
		}
	private: System::Windows::Forms::GroupBox^  ImageGroupBox;

	private: System::Windows::Forms::PictureBox^  ImageBox;

	private: System::ComponentModel::IContainer^  components;

	private:
		/// <summary>
		/// Required designer variable.
		/// </summary>


#pragma region Windows Form Designer generated code
		/// <summary>
		/// Required method for Designer support - do not modify
		/// the contents of this method with the code editor.
		/// </summary>
		void InitializeComponent(void)
		{
			this->components = (gcnew System::ComponentModel::Container());
			this->ImageGroupBox = (gcnew System::Windows::Forms::GroupBox());
			this->YMidLabel = (gcnew System::Windows::Forms::Label());
			this->XMidLabel = (gcnew System::Windows::Forms::Label());
			this->YMaxLabel = (gcnew System::Windows::Forms::Label());
			this->XMinLabel = (gcnew System::Windows::Forms::Label());
			this->YMinLabel = (gcnew System::Windows::Forms::Label());
			this->XMaxLabel = (gcnew System::Windows::Forms::Label());
			this->ImageBox = (gcnew System::Windows::Forms::PictureBox());
			this->starterTimer = (gcnew System::Windows::Forms::Timer(this->components));
			this->toolTip = (gcnew System::Windows::Forms::ToolTip(this->components));
			this->label13 = (gcnew System::Windows::Forms::Label());
			this->stepsBox = (gcnew System::Windows::Forms::NumericUpDown());
			this->alphaInputBox = (gcnew System::Windows::Forms::NumericUpDown());
			this->label1 = (gcnew System::Windows::Forms::Label());
			this->MInputBox = (gcnew System::Windows::Forms::NumericUpDown());
			this->Gamma0InputBox = (gcnew System::Windows::Forms::NumericUpDown());
			this->label2 = (gcnew System::Windows::Forms::Label());
			this->label3 = (gcnew System::Windows::Forms::Label());
			this->drawPhiButton = (gcnew System::Windows::Forms::Button());
			this->drawPsiButton = (gcnew System::Windows::Forms::Button());
			this->drawCpButton = (gcnew System::Windows::Forms::Button());
			this->flowOnlyButton = (gcnew System::Windows::Forms::Button());
			this->label5 = (gcnew System::Windows::Forms::Label());
			this->label6 = (gcnew System::Windows::Forms::Label());
			this->label7 = (gcnew System::Windows::Forms::Label());
			this->label8 = (gcnew System::Windows::Forms::Label());
			this->label9 = (gcnew System::Windows::Forms::Label());
			this->label10 = (gcnew System::Windows::Forms::Label());
			this->phiLabel = (gcnew System::Windows::Forms::Label());
			this->psiLabel = (gcnew System::Windows::Forms::Label());
			this->VLabel = (gcnew System::Windows::Forms::Label());
			this->CpLabel = (gcnew System::Windows::Forms::Label());
			this->xLabel = (gcnew System::Windows::Forms::Label());
			this->yLabel = (gcnew System::Windows::Forms::Label());
			this->ControlsGroupBox = (gcnew System::Windows::Forms::GroupBox());
			this->AutoStepsButton = (gcnew System::Windows::Forms::Button());
			this->RestartButton = (gcnew System::Windows::Forms::Button());
			this->NextStepButton = (gcnew System::Windows::Forms::Button());
			this->timeLabel = (gcnew System::Windows::Forms::Label());
			this->label11 = (gcnew System::Windows::Forms::Label());
			this->VScalButton = (gcnew System::Windows::Forms::Button());
			this->VScalLabel = (gcnew System::Windows::Forms::Label());
			this->label12 = (gcnew System::Windows::Forms::Label());
			this->animationTimer = (gcnew System::Windows::Forms::Timer(this->components));
			this->ImageGroupBox->SuspendLayout();
			(cli::safe_cast<System::ComponentModel::ISupportInitialize^>(this->ImageBox))->BeginInit();
			(cli::safe_cast<System::ComponentModel::ISupportInitialize^>(this->stepsBox))->BeginInit();
			(cli::safe_cast<System::ComponentModel::ISupportInitialize^>(this->alphaInputBox))->BeginInit();
			(cli::safe_cast<System::ComponentModel::ISupportInitialize^>(this->MInputBox))->BeginInit();
			(cli::safe_cast<System::ComponentModel::ISupportInitialize^>(this->Gamma0InputBox))->BeginInit();
			this->ControlsGroupBox->SuspendLayout();
			this->SuspendLayout();
			// 
			// ImageGroupBox
			// 
			this->ImageGroupBox->Anchor = static_cast<System::Windows::Forms::AnchorStyles>((((System::Windows::Forms::AnchorStyles::Top | System::Windows::Forms::AnchorStyles::Bottom)
				| System::Windows::Forms::AnchorStyles::Left)
				| System::Windows::Forms::AnchorStyles::Right));
			this->ImageGroupBox->Controls->Add(this->YMidLabel);
			this->ImageGroupBox->Controls->Add(this->XMidLabel);
			this->ImageGroupBox->Controls->Add(this->YMaxLabel);
			this->ImageGroupBox->Controls->Add(this->XMinLabel);
			this->ImageGroupBox->Controls->Add(this->YMinLabel);
			this->ImageGroupBox->Controls->Add(this->XMaxLabel);
			this->ImageGroupBox->Controls->Add(this->ImageBox);
			this->ImageGroupBox->Location = System::Drawing::Point(203, 6);
			this->ImageGroupBox->Name = L"ImageGroupBox";
			this->ImageGroupBox->Size = System::Drawing::Size(978, 658);
			this->ImageGroupBox->TabIndex = 0;
			this->ImageGroupBox->TabStop = false;
			// 
			// YMidLabel
			// 
			this->YMidLabel->AutoSize = true;
			this->YMidLabel->Location = System::Drawing::Point(34, 312);
			this->YMidLabel->Name = L"YMidLabel";
			this->YMidLabel->Size = System::Drawing::Size(10, 13);
			this->YMidLabel->TabIndex = 7;
			this->YMidLabel->Text = L"-";
			// 
			// XMidLabel
			// 
			this->XMidLabel->AutoSize = true;
			this->XMidLabel->Location = System::Drawing::Point(474, 622);
			this->XMidLabel->Name = L"XMidLabel";
			this->XMidLabel->Size = System::Drawing::Size(10, 13);
			this->XMidLabel->TabIndex = 6;
			this->XMidLabel->Text = L"-";
			// 
			// YMaxLabel
			// 
			this->YMaxLabel->AutoSize = true;
			this->YMaxLabel->Location = System::Drawing::Point(34, 16);
			this->YMaxLabel->Name = L"YMaxLabel";
			this->YMaxLabel->Size = System::Drawing::Size(10, 13);
			this->YMaxLabel->TabIndex = 5;
			this->YMaxLabel->Text = L"-";
			// 
			// XMinLabel
			// 
			this->XMinLabel->Anchor = static_cast<System::Windows::Forms::AnchorStyles>((System::Windows::Forms::AnchorStyles::Bottom | System::Windows::Forms::AnchorStyles::Left));
			this->XMinLabel->AutoSize = true;
			this->XMinLabel->Location = System::Drawing::Point(47, 622);
			this->XMinLabel->Name = L"XMinLabel";
			this->XMinLabel->Size = System::Drawing::Size(10, 13);
			this->XMinLabel->TabIndex = 4;
			this->XMinLabel->Text = L"-";
			// 
			// YMinLabel
			// 
			this->YMinLabel->Anchor = static_cast<System::Windows::Forms::AnchorStyles>((System::Windows::Forms::AnchorStyles::Bottom | System::Windows::Forms::AnchorStyles::Left));
			this->YMinLabel->AutoSize = true;
			this->YMinLabel->Location = System::Drawing::Point(34, 610);
			this->YMinLabel->Name = L"YMinLabel";
			this->YMinLabel->Size = System::Drawing::Size(10, 13);
			this->YMinLabel->TabIndex = 3;
			this->YMinLabel->Text = L"-";
			// 
			// XMaxLabel
			// 
			this->XMaxLabel->Anchor = static_cast<System::Windows::Forms::AnchorStyles>((System::Windows::Forms::AnchorStyles::Bottom | System::Windows::Forms::AnchorStyles::Right));
			this->XMaxLabel->AutoSize = true;
			this->XMaxLabel->Location = System::Drawing::Point(919, 622);
			this->XMaxLabel->Name = L"XMaxLabel";
			this->XMaxLabel->Size = System::Drawing::Size(10, 13);
			this->XMaxLabel->TabIndex = 2;
			this->XMaxLabel->Text = L"-";
			// 
			// ImageBox
			// 
			this->ImageBox->Anchor = static_cast<System::Windows::Forms::AnchorStyles>((((System::Windows::Forms::AnchorStyles::Top | System::Windows::Forms::AnchorStyles::Bottom)
				| System::Windows::Forms::AnchorStyles::Left)
				| System::Windows::Forms::AnchorStyles::Right));
			this->ImageBox->BackColor = System::Drawing::Color::White;
			this->ImageBox->BorderStyle = System::Windows::Forms::BorderStyle::FixedSingle;
			this->ImageBox->Cursor = System::Windows::Forms::Cursors::Cross;
			this->ImageBox->Location = System::Drawing::Point(50, 19);
			this->ImageBox->Name = L"ImageBox";
			this->ImageBox->Size = System::Drawing::Size(859, 600);
			this->ImageBox->TabIndex = 1;
			this->ImageBox->TabStop = false;
			this->ImageBox->SizeChanged += gcnew System::EventHandler(this, &MainForm::ImageBox_SizeChanged);
			this->ImageBox->Click += gcnew System::EventHandler(this, &MainForm::ImageBox_Click);
			this->ImageBox->MouseDown += gcnew System::Windows::Forms::MouseEventHandler(this, &MainForm::ImageBox_MouseDown);
			this->ImageBox->MouseEnter += gcnew System::EventHandler(this, &MainForm::ImageBox_MouseEnter);
			this->ImageBox->MouseLeave += gcnew System::EventHandler(this, &MainForm::ImageBox_MouseLeave);
			this->ImageBox->MouseMove += gcnew System::Windows::Forms::MouseEventHandler(this, &MainForm::ImageBox_MouseMove);
			this->ImageBox->MouseUp += gcnew System::Windows::Forms::MouseEventHandler(this, &MainForm::ImageBox_MouseUp);
			this->ImageBox->MouseWheel += gcnew System::Windows::Forms::MouseEventHandler(this, &MainForm::ImageBox_MouseWheel);
			// 
			// starterTimer
			// 
			this->starterTimer->Enabled = true;
			this->starterTimer->Interval = 1000;
			this->starterTimer->Tick += gcnew System::EventHandler(this, &MainForm::redrawTimer_Tick);
			// 
			// label13
			// 
			this->label13->AutoSize = true;
			this->label13->Location = System::Drawing::Point(21, 145);
			this->label13->Name = L"label13";
			this->label13->Size = System::Drawing::Size(80, 13);
			this->label13->TabIndex = 37;
			this->label13->Text = L"Автом. кроків:";
			this->toolTip->SetToolTip(this->label13, L"Кількість кроків анімації");
			// 
			// stepsBox
			// 
			this->stepsBox->Increment = System::Decimal(gcnew cli::array< System::Int32 >(4) { 10, 0, 0, 0 });
			this->stepsBox->Location = System::Drawing::Point(107, 143);
			this->stepsBox->Maximum = System::Decimal(gcnew cli::array< System::Int32 >(4) { 10000, 0, 0, 0 });
			this->stepsBox->Name = L"stepsBox";
			this->stepsBox->Size = System::Drawing::Size(87, 20);
			this->stepsBox->TabIndex = 36;
			this->toolTip->SetToolTip(this->stepsBox, L"Кількість кроків анімації");
			this->stepsBox->Value = System::Decimal(gcnew cli::array< System::Int32 >(4) { 1000, 0, 0, 0 });
			// 
			// alphaInputBox
			// 
			this->alphaInputBox->Increment = System::Decimal(gcnew cli::array< System::Int32 >(4) { 10, 0, 0, 0 });
			this->alphaInputBox->Location = System::Drawing::Point(107, 20);
			this->alphaInputBox->Maximum = System::Decimal(gcnew cli::array< System::Int32 >(4) { 180, 0, 0, 0 });
			this->alphaInputBox->Minimum = System::Decimal(gcnew cli::array< System::Int32 >(4) { 180, 0, 0, System::Int32::MinValue });
			this->alphaInputBox->Name = L"alphaInputBox";
			this->alphaInputBox->Size = System::Drawing::Size(80, 20);
			this->alphaInputBox->TabIndex = 0;
			this->alphaInputBox->ValueChanged += gcnew System::EventHandler(this, &MainForm::alphaInputBox_ValueChanged);
			// 
			// label1
			// 
			this->label1->AutoSize = true;
			this->label1->Location = System::Drawing::Point(44, 22);
			this->label1->Name = L"label1";
			this->label1->Size = System::Drawing::Size(24, 13);
			this->label1->TabIndex = 1;
			this->label1->Text = L"Кут";
			// 
			// MInputBox
			// 
			this->MInputBox->Increment = System::Decimal(gcnew cli::array< System::Int32 >(4) { 5, 0, 0, 0 });
			this->MInputBox->Location = System::Drawing::Point(107, 59);
			this->MInputBox->Minimum = System::Decimal(gcnew cli::array< System::Int32 >(4) { 4, 0, 0, 0 });
			this->MInputBox->Name = L"MInputBox";
			this->MInputBox->Size = System::Drawing::Size(80, 20);
			this->MInputBox->TabIndex = 2;
			this->MInputBox->Value = System::Decimal(gcnew cli::array< System::Int32 >(4) { 60, 0, 0, 0 });
			this->MInputBox->ValueChanged += gcnew System::EventHandler(this, &MainForm::MInputBox_ValueChanged);
			// 
			// Gamma0InputBox
			// 
			this->Gamma0InputBox->DecimalPlaces = 2;
			this->Gamma0InputBox->Enabled = false;
			this->Gamma0InputBox->Increment = System::Decimal(gcnew cli::array< System::Int32 >(4) { 1, 0, 0, 65536 });
			this->Gamma0InputBox->Location = System::Drawing::Point(107, 98);
			this->Gamma0InputBox->Minimum = System::Decimal(gcnew cli::array< System::Int32 >(4) { 100, 0, 0, System::Int32::MinValue });
			this->Gamma0InputBox->Name = L"Gamma0InputBox";
			this->Gamma0InputBox->Size = System::Drawing::Size(80, 20);
			this->Gamma0InputBox->TabIndex = 3;
			this->Gamma0InputBox->Value = System::Decimal(gcnew cli::array< System::Int32 >(4) { 1, 0, 0, 0 });
			this->Gamma0InputBox->ValueChanged += gcnew System::EventHandler(this, &MainForm::Gamma0InputBox_ValueChanged);
			// 
			// label2
			// 
			this->label2->AutoSize = true;
			this->label2->Location = System::Drawing::Point(15, 61);
			this->label2->Name = L"label2";
			this->label2->Size = System::Drawing::Size(84, 13);
			this->label2->TabIndex = 4;
			this->label2->Text = L"Кількість точок";
			// 
			// label3
			// 
			this->label3->AutoSize = true;
			this->label3->Location = System::Drawing::Point(25, 100);
			this->label3->Name = L"label3";
			this->label3->Size = System::Drawing::Size(64, 13);
			this->label3->TabIndex = 5;
			this->label3->Text = L"Циркуляція";
			// 
			// drawPhiButton
			// 
			this->drawPhiButton->Location = System::Drawing::Point(6, 350);
			this->drawPhiButton->Name = L"drawPhiButton";
			this->drawPhiButton->Size = System::Drawing::Size(182, 37);
			this->drawPhiButton->TabIndex = 9;
			this->drawPhiButton->Text = L"Потенціал течії";
			this->drawPhiButton->UseVisualStyleBackColor = true;
			this->drawPhiButton->Click += gcnew System::EventHandler(this, &MainForm::testButton_Click);
			// 
			// drawPsiButton
			// 
			this->drawPsiButton->Location = System::Drawing::Point(6, 392);
			this->drawPsiButton->Name = L"drawPsiButton";
			this->drawPsiButton->Size = System::Drawing::Size(182, 37);
			this->drawPsiButton->TabIndex = 10;
			this->drawPsiButton->Text = L"Лінії течії";
			this->drawPsiButton->UseVisualStyleBackColor = true;
			this->drawPsiButton->Click += gcnew System::EventHandler(this, &MainForm::drawPsiButton_Click);
			// 
			// drawCpButton
			// 
			this->drawCpButton->Location = System::Drawing::Point(6, 435);
			this->drawCpButton->Name = L"drawCpButton";
			this->drawCpButton->Size = System::Drawing::Size(182, 37);
			this->drawCpButton->TabIndex = 11;
			this->drawCpButton->Text = L"Тиск";
			this->drawCpButton->UseVisualStyleBackColor = true;
			this->drawCpButton->Click += gcnew System::EventHandler(this, &MainForm::drawCpButton_Click);
			// 
			// flowOnlyButton
			// 
			this->flowOnlyButton->Location = System::Drawing::Point(6, 517);
			this->flowOnlyButton->Name = L"flowOnlyButton";
			this->flowOnlyButton->Size = System::Drawing::Size(182, 37);
			this->flowOnlyButton->TabIndex = 12;
			this->flowOnlyButton->Text = L"Векторне поле течії";
			this->flowOnlyButton->UseVisualStyleBackColor = true;
			this->flowOnlyButton->Click += gcnew System::EventHandler(this, &MainForm::flowOnlyButton_Click);
			// 
			// label5
			// 
			this->label5->AutoSize = true;
			this->label5->Location = System::Drawing::Point(85, 166);
			this->label5->Name = L"label5";
			this->label5->Size = System::Drawing::Size(17, 13);
			this->label5->TabIndex = 14;
			this->label5->Text = L"X:";
			this->label5->Click += gcnew System::EventHandler(this, &MainForm::label5_Click);
			// 
			// label6
			// 
			this->label6->AutoSize = true;
			this->label6->Location = System::Drawing::Point(85, 181);
			this->label6->Name = L"label6";
			this->label6->Size = System::Drawing::Size(17, 13);
			this->label6->TabIndex = 15;
			this->label6->Text = L"Y:";
			this->label6->Click += gcnew System::EventHandler(this, &MainForm::label6_Click);
			// 
			// label7
			// 
			this->label7->AutoSize = true;
			this->label7->Location = System::Drawing::Point(41, 195);
			this->label7->Name = L"label7";
			this->label7->Size = System::Drawing::Size(61, 13);
			this->label7->TabIndex = 16;
			this->label7->Text = L"Потенціал:";
			this->label7->Click += gcnew System::EventHandler(this, &MainForm::label7_Click);
			// 
			// label8
			// 
			this->label8->AutoSize = true;
			this->label8->Location = System::Drawing::Point(26, 222);
			this->label8->Name = L"label8";
			this->label8->Size = System::Drawing::Size(76, 13);
			this->label8->TabIndex = 17;
			this->label8->Text = L"Функція течії:";
			this->label8->Click += gcnew System::EventHandler(this, &MainForm::label8_Click);
			// 
			// label9
			// 
			this->label9->AutoSize = true;
			this->label9->Location = System::Drawing::Point(2, 307);
			this->label9->Name = L"label9";
			this->label9->Size = System::Drawing::Size(100, 13);
			this->label9->TabIndex = 18;
			this->label9->Text = L"Вектор-швидкість:";
			this->label9->Click += gcnew System::EventHandler(this, &MainForm::label9_Click);
			// 
			// label10
			// 
			this->label10->AutoSize = true;
			this->label10->Location = System::Drawing::Point(67, 249);
			this->label10->Name = L"label10";
			this->label10->Size = System::Drawing::Size(35, 13);
			this->label10->TabIndex = 19;
			this->label10->Text = L"Тиск:";
			this->label10->Click += gcnew System::EventHandler(this, &MainForm::label10_Click);
			// 
			// phiLabel
			// 
			this->phiLabel->AutoSize = true;
			this->phiLabel->Location = System::Drawing::Point(104, 195);
			this->phiLabel->Name = L"phiLabel";
			this->phiLabel->Size = System::Drawing::Size(10, 13);
			this->phiLabel->TabIndex = 20;
			this->phiLabel->Text = L"-";
			this->phiLabel->Click += gcnew System::EventHandler(this, &MainForm::phiLabel_Click);
			// 
			// psiLabel
			// 
			this->psiLabel->AutoSize = true;
			this->psiLabel->Location = System::Drawing::Point(104, 222);
			this->psiLabel->Name = L"psiLabel";
			this->psiLabel->Size = System::Drawing::Size(10, 13);
			this->psiLabel->TabIndex = 21;
			this->psiLabel->Text = L"-";
			this->psiLabel->Click += gcnew System::EventHandler(this, &MainForm::psiLabel_Click);
			// 
			// VLabel
			// 
			this->VLabel->AutoSize = true;
			this->VLabel->Location = System::Drawing::Point(104, 307);
			this->VLabel->Name = L"VLabel";
			this->VLabel->Size = System::Drawing::Size(10, 13);
			this->VLabel->TabIndex = 22;
			this->VLabel->Text = L"-";
			this->VLabel->Click += gcnew System::EventHandler(this, &MainForm::VLabel_Click);
			// 
			// CpLabel
			// 
			this->CpLabel->AutoSize = true;
			this->CpLabel->Location = System::Drawing::Point(104, 249);
			this->CpLabel->Name = L"CpLabel";
			this->CpLabel->Size = System::Drawing::Size(10, 13);
			this->CpLabel->TabIndex = 23;
			this->CpLabel->Text = L"-";
			this->CpLabel->Click += gcnew System::EventHandler(this, &MainForm::CpLabel_Click);
			// 
			// xLabel
			// 
			this->xLabel->AutoSize = true;
			this->xLabel->Location = System::Drawing::Point(104, 166);
			this->xLabel->Name = L"xLabel";
			this->xLabel->Size = System::Drawing::Size(10, 13);
			this->xLabel->TabIndex = 24;
			this->xLabel->Text = L"-";
			this->xLabel->Click += gcnew System::EventHandler(this, &MainForm::xLabel_Click);
			// 
			// yLabel
			// 
			this->yLabel->AutoSize = true;
			this->yLabel->Location = System::Drawing::Point(104, 181);
			this->yLabel->Name = L"yLabel";
			this->yLabel->Size = System::Drawing::Size(10, 13);
			this->yLabel->TabIndex = 25;
			this->yLabel->Text = L"-";
			this->yLabel->Click += gcnew System::EventHandler(this, &MainForm::yLabel_Click);
			// 
			// ControlsGroupBox
			// 
			this->ControlsGroupBox->Anchor = static_cast<System::Windows::Forms::AnchorStyles>(((System::Windows::Forms::AnchorStyles::Top | System::Windows::Forms::AnchorStyles::Bottom)
				| System::Windows::Forms::AnchorStyles::Right));
			this->ControlsGroupBox->Controls->Add(this->AutoStepsButton);
			this->ControlsGroupBox->Controls->Add(this->RestartButton);
			this->ControlsGroupBox->Controls->Add(this->NextStepButton);
			this->ControlsGroupBox->Controls->Add(this->label13);
			this->ControlsGroupBox->Controls->Add(this->stepsBox);
			this->ControlsGroupBox->Controls->Add(this->timeLabel);
			this->ControlsGroupBox->Controls->Add(this->label11);
			this->ControlsGroupBox->Controls->Add(this->VScalButton);
			this->ControlsGroupBox->Controls->Add(this->VScalLabel);
			this->ControlsGroupBox->Controls->Add(this->label12);
			this->ControlsGroupBox->Controls->Add(this->yLabel);
			this->ControlsGroupBox->Controls->Add(this->xLabel);
			this->ControlsGroupBox->Controls->Add(this->CpLabel);
			this->ControlsGroupBox->Controls->Add(this->VLabel);
			this->ControlsGroupBox->Controls->Add(this->psiLabel);
			this->ControlsGroupBox->Controls->Add(this->phiLabel);
			this->ControlsGroupBox->Controls->Add(this->label10);
			this->ControlsGroupBox->Controls->Add(this->label9);
			this->ControlsGroupBox->Controls->Add(this->label8);
			this->ControlsGroupBox->Controls->Add(this->label7);
			this->ControlsGroupBox->Controls->Add(this->label6);
			this->ControlsGroupBox->Controls->Add(this->label5);
			this->ControlsGroupBox->Controls->Add(this->flowOnlyButton);
			this->ControlsGroupBox->Controls->Add(this->drawCpButton);
			this->ControlsGroupBox->Controls->Add(this->drawPsiButton);
			this->ControlsGroupBox->Controls->Add(this->drawPhiButton);
			this->ControlsGroupBox->Controls->Add(this->label3);
			this->ControlsGroupBox->Controls->Add(this->label2);
			this->ControlsGroupBox->Controls->Add(this->Gamma0InputBox);
			this->ControlsGroupBox->Controls->Add(this->MInputBox);
			this->ControlsGroupBox->Controls->Add(this->label1);
			this->ControlsGroupBox->Controls->Add(this->alphaInputBox);
			this->ControlsGroupBox->Location = System::Drawing::Point(12, 6);
			this->ControlsGroupBox->Name = L"ControlsGroupBox";
			this->ControlsGroupBox->Size = System::Drawing::Size(197, 658);
			this->ControlsGroupBox->TabIndex = 1;
			this->ControlsGroupBox->TabStop = false;
			this->ControlsGroupBox->Enter += gcnew System::EventHandler(this, &MainForm::ControlsGroupBox_Enter);
			// 
			// AutoStepsButton
			// 
			this->AutoStepsButton->Location = System::Drawing::Point(101, 559);
			this->AutoStepsButton->Name = L"AutoStepsButton";
			this->AutoStepsButton->Size = System::Drawing::Size(87, 38);
			this->AutoStepsButton->TabIndex = 40;
			this->AutoStepsButton->Text = L"Пройти кроки";
			this->AutoStepsButton->UseVisualStyleBackColor = true;
			this->AutoStepsButton->Click += gcnew System::EventHandler(this, &MainForm::AutoStepsButton_Click);
			// 
			// RestartButton
			// 
			this->RestartButton->Location = System::Drawing::Point(6, 559);
			this->RestartButton->Name = L"RestartButton";
			this->RestartButton->Size = System::Drawing::Size(89, 38);
			this->RestartButton->TabIndex = 39;
			this->RestartButton->Text = L"Заново";
			this->RestartButton->UseVisualStyleBackColor = true;
			this->RestartButton->Click += gcnew System::EventHandler(this, &MainForm::RestartButton_Click);
			// 
			// NextStepButton
			// 
			this->NextStepButton->Location = System::Drawing::Point(6, 603);
			this->NextStepButton->Name = L"NextStepButton";
			this->NextStepButton->Size = System::Drawing::Size(182, 38);
			this->NextStepButton->TabIndex = 38;
			this->NextStepButton->Text = L"Наступний крок >>>";
			this->NextStepButton->UseVisualStyleBackColor = true;
			this->NextStepButton->Click += gcnew System::EventHandler(this, &MainForm::NextStepButton_Click);
			// 
			// timeLabel
			// 
			this->timeLabel->AutoSize = true;
			this->timeLabel->Location = System::Drawing::Point(104, 121);
			this->timeLabel->Name = L"timeLabel";
			this->timeLabel->Size = System::Drawing::Size(22, 13);
			this->timeLabel->TabIndex = 34;
			this->timeLabel->Text = L"0.0";
			// 
			// label11
			// 
			this->label11->AutoSize = true;
			this->label11->Location = System::Drawing::Point(71, 121);
			this->label11->Name = L"label11";
			this->label11->Size = System::Drawing::Size(30, 13);
			this->label11->TabIndex = 33;
			this->label11->Text = L"Час:";
			// 
			// VScalButton
			// 
			this->VScalButton->Location = System::Drawing::Point(6, 477);
			this->VScalButton->Name = L"VScalButton";
			this->VScalButton->Size = System::Drawing::Size(182, 37);
			this->VScalButton->TabIndex = 29;
			this->VScalButton->Text = L"Швидкість";
			this->VScalButton->UseVisualStyleBackColor = true;
			this->VScalButton->Click += gcnew System::EventHandler(this, &MainForm::VScalButton_Click);
			// 
			// VScalLabel
			// 
			this->VScalLabel->AutoSize = true;
			this->VScalLabel->Location = System::Drawing::Point(104, 277);
			this->VScalLabel->Name = L"VScalLabel";
			this->VScalLabel->Size = System::Drawing::Size(10, 13);
			this->VScalLabel->TabIndex = 28;
			this->VScalLabel->Text = L"-";
			this->VScalLabel->Click += gcnew System::EventHandler(this, &MainForm::VScalLabel_Click);
			// 
			// label12
			// 
			this->label12->AutoSize = true;
			this->label12->Location = System::Drawing::Point(40, 277);
			this->label12->Name = L"label12";
			this->label12->Size = System::Drawing::Size(62, 13);
			this->label12->TabIndex = 27;
			this->label12->Text = L"Швидкість:";
			this->label12->Click += gcnew System::EventHandler(this, &MainForm::label12_Click);
			// 
			// animationTimer
			// 
			this->animationTimer->Interval = 200;
			this->animationTimer->Tick += gcnew System::EventHandler(this, &MainForm::animationTimer_Tick);
			// 
			// MainForm
			// 
			this->AutoScaleDimensions = System::Drawing::SizeF(6, 13);
			this->AutoScaleMode = System::Windows::Forms::AutoScaleMode::Font;
			this->ClientSize = System::Drawing::Size(1184, 661);
			this->Controls->Add(this->ControlsGroupBox);
			this->Controls->Add(this->ImageGroupBox);
			this->Name = L"MainForm";
			this->Load += gcnew System::EventHandler(this, &MainForm::MainForm_Load);
			this->ImageGroupBox->ResumeLayout(false);
			this->ImageGroupBox->PerformLayout();
			(cli::safe_cast<System::ComponentModel::ISupportInitialize^>(this->ImageBox))->EndInit();
			(cli::safe_cast<System::ComponentModel::ISupportInitialize^>(this->stepsBox))->EndInit();
			(cli::safe_cast<System::ComponentModel::ISupportInitialize^>(this->alphaInputBox))->EndInit();
			(cli::safe_cast<System::ComponentModel::ISupportInitialize^>(this->MInputBox))->EndInit();
			(cli::safe_cast<System::ComponentModel::ISupportInitialize^>(this->Gamma0InputBox))->EndInit();
			this->ControlsGroupBox->ResumeLayout(false);
			this->ControlsGroupBox->PerformLayout();
			this->ResumeLayout(false);

		}
#pragma endregion

	private: Vector2f pixelToReal(std::pair<int, int> pixel)
			 {
				 Vector2f res;

				 pixel.second = ImageBox->Height - pixel.second;

				 res.x = pixel.first * 1.0f / zoom;
				 res.y = pixel.second * 1.0f / zoom;
				 res = res + (*offset);

				 return res;
			 }

	private: float pixelToReal(int pixel)
			 {
				 return pixel * 1.0f / zoom;
			 }

	private: std::pair<int, int> realToPixel(Vector2f real)
			 {
				 std::pair<int, int> res;

				 real = real - (*offset);

				 res.first  = int(real.x * zoom);
				 res.second = int(real.y * zoom);

				 res.second = ImageBox->Height - res.second;

				 return res;
			 }

	private: int realToPixel(float real)
			 {
				 return int(real * zoom);
			 }

	private: void drawCurve()
			 {
				 Graphics^ g = bufferedGraphics->Graphics;
				 Pen^ pen = gcnew Pen(Color::Blue, 2);

				 std::pair<int, int> p2, p3, p4;
				 p2 = realToPixel(Vector2f(-2.0f / 9.0f, 2.0f / 9.0f));
				 p3 = realToPixel(Vector2f(0.0f, -2.0f / 9.0f));
				 p4 = realToPixel(Vector2f(2.0f / 9.0f, 2.0f / 9.0f));

				 g->DrawLine(pen, p2.first, p2.second, p3.first, p3.second);
				 g->DrawLine(pen, p3.first, p3.second, p4.first, p4.second);

			 
			 }

	private: void drawGrid()
			 {
				 Graphics^ g = bufferedGraphics->Graphics;
				 Color gridColor = Color::FromArgb(100, Color::LightBlue);

				 Pen^ pen = gcnew Pen(gridColor);

				 for (int i = 1; i < 10; i++)
					 g->DrawLine(pen,
					 i * ImageBox->Width / 10,
					 0,
					 i * ImageBox->Width / 10,
					 ImageBox->Height);

				 for (int i = 1; i < 10; i++)
					 g->DrawLine(pen,
					 0,
					 i * ImageBox->Height / 10,
					 ImageBox->Width,
					 i * ImageBox->Height / 10);
			 }

	private: void drawTails()
	{
		const float dotRadius = 0.002f;

		Graphics^ g = bufferedGraphics->Graphics;
		Pen^ pen = gcnew Pen(Color::Blue);
		SolidBrush^ brush = gcnew SolidBrush(Color::Blue);
		SolidBrush^ brush1 = gcnew SolidBrush(Color::Blue);//Green);

		std::vector< std::vector<Vector2f> > tails = solver->getTails();
		std::vector<Vector2f> tailSources = solver->getTailSources();
		std::vector< std::vector<float> > tailGamma = solver->getTailGamma();

		std::pair<int, int> p1, p2;

		/*for (int tail = 0; tail < int(tails.size()); tail++)
		{
			for (int i = 0; i < int(tails[tail].size()) - 1; i++)
			{
				p1 = realToPixel(tails[tail][i]);
				p2 = realToPixel(tails[tail][i + 1]);
				g->DrawLine(pen, p1.first, p1.second, p2.first, p2.second);
			}
			if (int(tails[tail].size()) > 0)
			{
				p1 = realToPixel(tails[tail][int(tails[tail].size()) - 1]);
				p2 = realToPixel(tailSources[tail]);
				g->DrawLine(pen, p1.first, p1.second, p2.first, p2.second);
			}
		}*/

		for (int tail = 0; tail < int(tails.size()); tail++)
			for (int i = 0; i < int(tails[tail].size()); i++)
			{
				if (tailGamma[tail][i] == 0)
					continue;
				p1 = realToPixel(tails[tail][i] + Vector2f(-dotRadius, dotRadius));
				if (tailGamma[tail][i] > 0)
				{
					if (i == 0)
						g->FillRectangle(brush, p1.first, p1.second,
							realToPixel(10 * dotRadius), realToPixel(10 * dotRadius));
					else
						g->FillRectangle(brush, p1.first, p1.second,
						realToPixel(5 * dotRadius), realToPixel(5 * dotRadius));
				}
				else
				{
					if (i == 0)
						g->FillRectangle(brush1, p1.first, p1.second,
						realToPixel(10 * dotRadius), realToPixel(10 * dotRadius));
					else
						g->FillRectangle(brush1, p1.first, p1.second,
						realToPixel(5 * dotRadius), realToPixel(5 * dotRadius));
				}
			}
	}

	private: void drawPoints()
			 {
				 /*const float dotRadius = 0.01f;
				 const int zoneAlpha = 100;
				 Color dotColor = Color::Crimson;

				 Graphics^ g = bufferedGraphics->Graphics;
				 SolidBrush^ brush = gcnew SolidBrush(dotColor);
				 Pen^ pen = gcnew Pen(dotColor);
				 SolidBrush^ zoneBrush = gcnew SolidBrush(
					 Color::FromArgb(zoneAlpha, dotColor));

				 float delta = solver->getDelta();
				 std::vector<Vector2f> singularities = solver->getSingularities();
				 std::vector<Vector2f> colocations = solver->getColocations();

				 for (int i = 0; i < int(singularities.size()); ++i)
				 {
					 std::pair<int, int> p = realToPixel(singularities[i] + 
						 Vector2f(-dotRadius, dotRadius));
					 g->FillEllipse(brush, p.first, p.second, 
						 realToPixel(2 * dotRadius), realToPixel(2 * dotRadius));
				 }

				 for (int i = 0; i < int(colocations.size()); ++i)
				 {
					 std::pair<int, int> p1, p2, p3, p4;
					 p1 = realToPixel(Vector2f(colocations[i].x - dotRadius / 2.0f, 
						 colocations[i].y - dotRadius / 2.0f));
					 p2 = realToPixel(Vector2f(colocations[i].x + dotRadius / 2.0f, 
						 colocations[i].y + dotRadius / 2.0f));
					 p3 = realToPixel(Vector2f(colocations[i].x - dotRadius / 2.0f, 
						 colocations[i].y + dotRadius / 2.0f));
					 p4 = realToPixel(Vector2f(colocations[i].x + dotRadius / 2.0f, 
						 colocations[i].y - dotRadius / 2.0f));

					 g->DrawLine(pen, p1.first, p1.second, p2.first, p2.second);
					 g->DrawLine(pen, p3.first, p3.second, p4.first, p4.second);
				 }*/
			 }

	private: void drawGradient(bool flag)
			 {
				 int gradientWidth = 80;
				 int gradientHeight = 560;
				 const int clasters = 24;
				 Color lowColor = Color::White;
				 Color highColor = Color::Black;

				 if (flag)
				 {
					 Graphics^ g = bufferedGraphics->Graphics;
					 Pen^ pen = gcnew Pen(Color::Black);
					 SolidBrush^ brush = gcnew SolidBrush(this->BackColor);

					 g->FillRectangle(brush, ImageBox->Width - gradientWidth - 5, 0,
						 gradientWidth + 5, gradientHeight + 5);

					 g->DrawRectangle(pen, ImageBox->Width - gradientWidth - 5, 0,
						 gradientWidth + 5, gradientHeight + 5);

					 for (int i = 0; i < clasters; i++)
					 {
						 int claster = clasters - i - 1;
						 float alpha = claster * 1.0f / clasters;
						 Color color = Color::FromArgb(
							 int(alpha * highColor.R + (1.0f - alpha) * lowColor.R),
							 int(alpha * highColor.G + (1.0f - alpha) * lowColor.G),
							 int(alpha * highColor.B + (1.0f - alpha) * lowColor.B));

						 brush->Color = color;

						 g->FillRectangle(brush, 
							 ImageBox->Width - gradientWidth, gradientHeight * i / clasters,
							 gradientWidth, gradientHeight / clasters);

						 g->DrawRectangle(pen, 
							 ImageBox->Width - gradientWidth, gradientHeight * i / clasters,
							 gradientWidth, gradientHeight / clasters);
					 }

					 for (int i = 0; i < clasters + 1; i++)
					 {
						 GradientLabels[i]->Text = floatToString(minValue + (maxValue-minValue) * (clasters - i) / clasters);
					 }
				 }
				 else
				 {
					 for (int i = 0; i < clasters + 1; i++)
					 {
						 GradientLabels[i]->Text = "";
					 }
				 }
			 }

	private: void drawScalarField(float (FlowSolver::*func)(float, float) const)
			 {
				 const int clasters = 24;
				 Color lowColor = Color::White;
				 Color highColor = Color::Black;

				 Graphics^ g = bufferedGraphics->Graphics;
				 Pen^ pen = gcnew Pen(Color::Black);

				 value->clear();
				 value->resize( ImageBox->Height, std::vector<float>(ImageBox->Width, 0.0f));
				 minValue = Inf;
				 maxValue = - Inf;

				 for (int x = 0; x < ImageBox->Width; x++)
				 {
					 for (int y = 0; y < ImageBox->Height; y++)
					 {
						 Vector2f p = pixelToReal(std::make_pair(x, y));
						 (*value)[y][x] = (solver->*func)(p.x, p.y);
						 float v = (*value)[y][x];
						 if (minValue > v)
							 minValue = v;
						 if (maxValue < v)
						 {
							 maxValue = v;
							 /*this->Text = Convert::ToString(p.x) + " " + Convert::ToString(p.y) +
							 " " + Convert::ToString(v);*/
						 }
					 }
				 }

				 for (int x = 0; x < ImageBox->Width; x++)
				 {
					 for (int y = 0; y < ImageBox->Height; y++)
					 {						 
						 int claster;
						 claster = int(clasters * ((*value)[y][x] - minValue) / (maxValue - minValue));
						 claster = std::max(std::min(claster, clasters - 1), 0);

						 float alpha = claster * 1.0f / clasters;
						 Color color = Color::FromArgb(
							 int(alpha * highColor.R + (1.0f - alpha) * lowColor.R),
							 int(alpha * highColor.G + (1.0f - alpha) * lowColor.G),
							 int(alpha * highColor.B + (1.0f - alpha) * lowColor.B));

						 pen->Color = color;
						 g->DrawLine(pen, x, y, x + 1, y + 1);
					 }
				 }	 
			 }

	private: void drawVectorField()
			 {
				 const int step = 10;
				 Color penColor = Color::Red;

				 Graphics^ g = bufferedGraphics->Graphics;
				 Pen^ pen = gcnew Pen(penColor);

				 for (int x = step / 2; x < ImageBox->Width; x += step)
				 {
					 for (int y = step / 2; y < ImageBox->Height; y += step)
					 {
						 Vector2f p = pixelToReal(std::make_pair(x, y));
						 Vector2f v = solver->V(p.x, p.y);

						 v = v * step;

						 g->DrawLine(pen, x, y, x + (int(v.x)), y - (int(v.y)));
					 }
				 }				 
			 }

	private: void redraw(int paint)
			 {
				 Graphics^ g = bufferedGraphics->Graphics;
				 g->Clear(Color::White);

				 if (paint == 6)
				 {
					 drawScalarField(&FlowSolver::phiVihr);
				 }
				 if (paint == 1)
				 {
					 drawScalarField(&FlowSolver::phiContinuous);
				 }
				 else if (paint == 2)
				 {
					 drawScalarField(&FlowSolver::psi);
				 }
				 else if (paint == 3)
				 {
					 drawScalarField(&FlowSolver::Cp);
				 }
				 else if (paint == 4)
				 {
					 drawScalarField(&FlowSolver::V_scal);
				 }

				 if (0 < paint && paint <= 6)
					 drawVectorField();

				 drawCurve();
				 drawPoints();
				 drawTails();
				 drawGrid();

				 if (0 < paint && paint <= 4 || paint == 6 )
					 drawGradient(true);
				 else
					 drawGradient(false);

				 graphicsState = g->Save();

				 bufferedGraphics->Render();
				 return;
			 }
	private: System::Void alphaInputBox_ValueChanged(System::Object^  sender, System::EventArgs^  e) {
				 alpha = float(Convert::ToDouble(alphaInputBox->Value));
				 alpha = alpha * Pi / 180.0f;
				 solver->setAlpha(alpha);
				 solver->init();

				 redraw(0);
			 }
	private: System::Void MInputBox_ValueChanged(System::Object^  sender, System::EventArgs^  e) {
				 M = Convert::ToInt32(MInputBox->Value);
				 solver->setM(M);
				 solver->init();

				 redraw(0);
			 }
	private: System::Void Gamma0InputBox_ValueChanged(System::Object^  sender, System::EventArgs^  e) {
				 gamma0 = float(Convert::ToDouble(Gamma0InputBox->Value));
				 solver->setGamma0(gamma0);
				 solver->init();

				 redraw(0);
			 }
	private: System::Void ImageBox_SizeChanged(System::Object^  sender, System::EventArgs^  e) {
				 if (ImageBox->Width == 0 || ImageBox->Height == 0)
					 return;

				 bufferedGraphicsContext->MaximumBuffer = Drawing::Size(ImageBox->Width, ImageBox->Height);
				 bufferedGraphics = bufferedGraphicsContext->Allocate(ImageBox->CreateGraphics(), 
					 System::Drawing::Rectangle(Point(0, 0), Drawing::Size(ImageBox->Width, ImageBox->Height)));

				 XMidLabel->Top = ImageBox->Bottom + 1;
				 XMidLabel->Left = (ImageBox->Left + ImageBox->Right) / 2 - XMidLabel->Width / 2 + 2;
				 YMidLabel->Top = ImageBox->Top + ImageBox->Height / 2 - 5;

				 redraw(0);
			 }
	private: void updateZoom(float zoom)
			 {
				 this->zoom = zoom;
			 }
	private: System::Void ImageBox_MouseWheel(System::Object^  sender, System::Windows::Forms::MouseEventArgs^  e) {
				 float new_zoom = float(int(zoom + float(e->Delta / 10) + 0.1));
				 if (int(new_zoom) < 10)
					 new_zoom = 10.0f;
				 if (int(new_zoom) > 5000)
					 new_zoom = 5000.0f;

				 offset->x += e->X * 1.0f / zoom;
				 offset->y += (ImageBox->Height - e->Y) * 1.0f / zoom;
				 offset->x -= e->X * 1.0f / new_zoom;
				 offset->y -= (ImageBox->Height - e->Y) * 1.0f / new_zoom;

				 updateZoom(new_zoom);

				 updateAxisLabels();
				 redraw(0);
			 }
	private: System::Void ImageBox_MouseEnter(System::Object^  sender, System::EventArgs^  e) {
				 ImageBox->Focus();
			 }
	private: System::Void redrawTimer_Tick(System::Object^  sender, System::EventArgs^  e) {
				 Graphics^ g = bufferedGraphics->Graphics;
				 g->Restore(graphicsState);
				 bufferedGraphics->Render();
			 }
	private: System::Void ImageBox_MouseMove(System::Object^  sender, System::Windows::Forms::MouseEventArgs^  e) {
				 if (isDragged)
				 {
					 offset->x -= (e->X - prevX) / zoom;
					 offset->y += (e->Y - prevY) / zoom;

					 prevX = e->X;
					 prevY = e->Y;

					 updateAxisLabels();
					 redraw(0);
				 }
				 else
				 {
					 Vector2f p = pixelToReal(std::make_pair(e->X, e->Y));

					 xLabel->Text = Convert::ToString(p.x);
					 yLabel->Text = Convert::ToString(p.y);

					 phiLabel->Text = Convert::ToString(solver->phi(p.x, p.y));
					 psiLabel->Text = Convert::ToString(solver->psi(p.x, p.y));
					 CpLabel->Text = Convert::ToString(solver->Cp(p.x, p.y));

					 Vector2f t = solver->V(p.x, p.y);
					 VScalLabel->Text = Convert::ToString(solver->V_scal(p.x, p.y));
					 VLabel->Text = "(" + Convert::ToString(t.x) + ";\n" + Convert::ToString(t.y) + ")";
				 }
			 }
	private: System::Void ImageBox_MouseDown(System::Object^  sender, System::Windows::Forms::MouseEventArgs^  e) {
				 isDragged = true;
				 ImageBox->Cursor = Cursors::SizeAll;
				 prevX = e->X;
				 prevY = e->Y;
			 }
	private: System::Void ImageBox_MouseUp(System::Object^  sender, System::Windows::Forms::MouseEventArgs^  e) {
				 isDragged = false;
				 ImageBox->Cursor = Cursors::Cross;
			 }
	private: System::Void drawPhiButtonVihr_Click(System::Object^  sender, System::EventArgs^  e) {
				 animationTimer->Enabled = false;
				 redraw(6);
			 }
	private: System::Void testButton_Click(System::Object^  sender, System::EventArgs^  e) {
			 	 animationTimer->Enabled = false;
				 redraw(1);
			 }
	private: System::Void drawPsiButton_Click(System::Object^  sender, System::EventArgs^  e) {
				 animationTimer->Enabled = false;
				 redraw(2);
			 }
	private: System::Void drawCpButton_Click(System::Object^  sender, System::EventArgs^  e) {
				 animationTimer->Enabled = false;
				 redraw(3);
			 }
	private: System::Void VScalButton_Click(System::Object^  sender, System::EventArgs^  e) {
				 animationTimer->Enabled = false;
				 redraw(4);
			 }
	private: System::Void flowOnlyButton_Click(System::Object^  sender, System::EventArgs^  e) {
				 animationTimer->Enabled = false;	
				 redraw(5);
			 }
	private: System::Void ImageBox_MouseLeave(System::Object^  sender, System::EventArgs^  e) {
				 xLabel->Text = "-";
				 yLabel->Text = "-";
				 phiLabel->Text = "-";
				 psiLabel->Text = "-";
				 VLabel->Text = "-";
				 VScalLabel->Text = "-";
				 CpLabel->Text = "-";
			 }
	private: String^ floatToString(float f)
			 {
				 String^ t = "";
				 int tpos;

				 t = Convert::ToString(f);
				 tpos = std::max(t->IndexOf(","), t->IndexOf("."));
				 if (tpos != -1 && t->Length > tpos + 4) 
					 t = t->Remove(tpos + 4);

				 return t;
			 }
	private: void updateAxisLabels() {
				 Vector2f LBCorner, RTCorner, Middle;

				 LBCorner = pixelToReal(std::make_pair(0, 0));
				 RTCorner = pixelToReal(std::make_pair(ImageBox->Width, ImageBox->Height));
				 Middle = pixelToReal(std::make_pair(ImageBox->Width / 2, ImageBox->Height / 2));

				 XMaxLabel->Text = floatToString(RTCorner.x);
				 XMaxLabel->Left = ImageBox->Right - XMaxLabel->Width / 2;

				 YMaxLabel->Text = floatToString(LBCorner.y);
				 YMaxLabel->Left = ImageBox->Left - YMaxLabel->Width;

				 XMinLabel->Text = floatToString(LBCorner.x);
				 XMinLabel->Left = ImageBox->Left - XMinLabel->Width / 2;

				 YMinLabel->Text = floatToString(RTCorner.y);
				 YMinLabel->Left = ImageBox->Left - YMinLabel->Width;

				 XMidLabel->Text = floatToString(Middle.x);
				 XMidLabel->Left = (ImageBox->Left + ImageBox->Right) / 2 - XMidLabel->Width / 2 + 2;

				 YMidLabel->Text = floatToString(Middle.y);
				 YMidLabel->Left = ImageBox->Left - YMidLabel->Width;
			 }
	private: System::Void redButton_Click(System::Object^  sender, System::EventArgs^  e) {
				 //float px = 0.0f;
				 //float py = 0.0f;
				 //float l, r;
				 //float v2, v3;
				 //l = -10.0f;
				 //r = +10.0f;
				 //while (r - l > 0.0001)
				 //{
					// float a, b;
					// a = l + (r - l) / 3.0f;
					// b = l + (r - l) / 3.0f * 2.0f;

					// Vector2f tmp;
					// solver->setGamma0(a); solver->solve(); tmp = solver->V(px, py); 
					// v2 = tmp.x * tmp.x + tmp.y * tmp.y;
					// solver->setGamma0(b); solver->solve(); tmp = solver->V(px, py); 
					// v3 = tmp.x * tmp.x + tmp.y * tmp.y;

					// if (v2 < v3)
					//	 r = b;
					// else
					//	 l = a;
				 //}

				 //solver->setGamma0(gamma0); solver->solve();
				 ///*MessageBox::Show("Gamma0 = " + Convert::ToString(l) + 
					// "\nV = " + Convert::ToString(sqrt(v2)));*/
				 //MessageBox::Show("Gamma0 = " + Convert::ToString(l));

				 float px = 0.0f;
				 float py = 0.5f;
				 float mg = 0.0f, mv = Inf;
				 float v;
				 for (float g = -10.0f; g <= 10.0f; g += 0.01f )
				 {
					 solver->setGamma0(g); 
					 solver->init(); 
					 v = (solver->V(px, py)).length();
					 if (v < mv)
					 {
						 mv = v;
						 mg = g;
					 }
				 }
				 MessageBox::Show("Gamma0 = " + Convert::ToString(mg) + 
					"\nV = " + Convert::ToString(sqrt(mv)));
			 }
	private: System::Void NextStepButton_Click(System::Object^  sender, System::EventArgs^  e) {
		animationTimer->Enabled = false;
		solver->nextStep();
		timeLabel->Text = floatToString(solver->getTime());
		redraw(0);
	}
	private: System::Void RestartButton_Click(System::Object^  sender, System::EventArgs^  e) {
		animationTimer->Enabled = false;

		solver->restart();

		timeLabel->Text = floatToString(solver->getTime());
		redraw(0);
	}
	private: System::Void animationTimer_Tick(System::Object^  sender, System::EventArgs^  e) {
		if (step >= Convert::ToInt32(stepsBox->Value))
		{
			animationTimer->Enabled = false;
			return;
		}
		else
		{
			step++;
			solver->nextStep();
			timeLabel->Text = floatToString(solver->getTime());
			redraw(5);
			//redraw(0);
		}
	}
	private: System::Void AutoStepsButton_Click(System::Object^  sender, System::EventArgs^  e) {
		step = 0;
		animationTimer->Enabled = true;
	}
	private: System::Void ImageBox_Click(System::Object^  sender, System::EventArgs^  e) {
			 }
private: System::Void MainForm_Load(System::Object^  sender, System::EventArgs^  e) {
		 }
private: System::Void ControlsGroupBox_Enter(System::Object^  sender, System::EventArgs^  e) {
		 }
private: System::Void label5_Click(System::Object^  sender, System::EventArgs^  e) {
		 }
private: System::Void xLabel_Click(System::Object^  sender, System::EventArgs^  e) {
		 }
private: System::Void yLabel_Click(System::Object^  sender, System::EventArgs^  e) {
		 }
private: System::Void label6_Click(System::Object^  sender, System::EventArgs^  e) {
		 }
private: System::Void label7_Click(System::Object^  sender, System::EventArgs^  e) {
		 }
private: System::Void phiLabel_Click(System::Object^  sender, System::EventArgs^  e) {
		 }
private: System::Void psiLabel_Click(System::Object^  sender, System::EventArgs^  e) {
		 }
private: System::Void label8_Click(System::Object^  sender, System::EventArgs^  e) {
		 }
private: System::Void label10_Click(System::Object^  sender, System::EventArgs^  e) {
		 }
private: System::Void CpLabel_Click(System::Object^  sender, System::EventArgs^  e) {
		 }
private: System::Void VScalLabel_Click(System::Object^  sender, System::EventArgs^  e) {
		 }
private: System::Void VLabel_Click(System::Object^  sender, System::EventArgs^  e) {
		 }
private: System::Void label9_Click(System::Object^  sender, System::EventArgs^  e) {
		 }
private: System::Void label12_Click(System::Object^  sender, System::EventArgs^  e) {
		 }
private: System::Void button1_Click(System::Object^  sender, System::EventArgs^  e) {
}
};
}
