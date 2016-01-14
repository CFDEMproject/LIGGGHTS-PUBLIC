#pragma once

#include "..\..\..\..\UtilityLib\tarray.h"
#include "..\..\..\..\UtilityLib\mystring.cpp"
#include "..\..\..\..\UtilityLib\myfile.cpp"
#include<vcclr.h>
#include <ctime>
#include "log_tool_classes.h"

namespace HOTINT_Log {

	using namespace System;
	using namespace System::ComponentModel;
	using namespace System::Collections;
	using namespace System::Windows::Forms;
	using namespace System::Data;
	using namespace System::Drawing;

#define IDENT_LOG true
#define IDENT_BUG false


//#define CHANGES_LOG_FILENAME "..\\..\\..\\..\\documentation\\changes_log.txt"
//#define CHANGES_LOG_FILENAME "..//..//..//..//documentation//changes_log.txt"
//#define CHANGES_LOG_FILENAME "../../../../documentation/changes_log.txt"
//#define CHANGES_LOG_FILENAME "../../../documentation/changes_log.txt"
//#define CHANGES_LOG_FILENAME "changes_log.txt"


//for final use
//#define USERNAME_FILENAME "../../userdata/LOG_USER.txt"
//#define LOG_FILENAME "../../userdata/LOG_Config.txt"	

//for test purposes, since the debug/release folder is 2 times deeper than the final use folder
//#define USERNAME_FILENAME "../../../../userdata/LOG_USER.txt"
//#define LOG_FILENAME "../../../../userdata/LOG_Config.txt"
#define USERNAME_FILENAME "../LOG_USER.txt"
#define LOG_FILENAME "../LOG_Config.txt"	
	


//#define USERNAME_FILENAME "LOG_USER.txt"
//#define USERNAME_FILENAME "LOG_Config.txt"


#define IDENT_START_LOG "CHANGES:"
#define IDENT_START_BUG "KNOWN/OPEN BUGS:"
#define IDENT_START_SOLVED_BUGS "SOLVED BUGS:"
#define IDENT_BUGS_ADD_SOLUTION "***PLEASE ADD DESCRIPTION OF SOLUTION***"
#define IDENT_SOLUTION "SOLUTION:"
#define IDENT_NEW_LOG "@@"
#define IDENT_NEW_BUG "@@B"


	/// <summary>
	/// Zusammenfassung für Form1
	///
	/// Warnung: Wenn Sie den Namen dieser Klasse ändern, müssen Sie auch
	///          die Ressourcendateiname-Eigenschaft für das Tool zur Kompilierung verwalteter Ressourcen ändern,
	///          das allen RESX-Dateien zugewiesen ist, von denen diese Klasse abhängt.
	///          Anderenfalls können die Designer nicht korrekt mit den lokalisierten Ressourcen
	///          arbeiten, die diesem Formular zugewiesen sind.
	/// </summary>
	public ref class Form1 : public System::Windows::Forms::Form
	{
	public:
		Form1(void)
		{
			InitializeComponent();
			//
			//TODO: Konstruktorcode hier hinzufügen.
			//
		}

	protected:
		/// <summary>
		/// Verwendete Ressourcen bereinigen.
		/// </summary>
		~Form1()
		{
			if (components)
			{
				delete components;
			}
			Destroy();
			//delete lastLogData;
			//delete header;
		}
	private: System::Windows::Forms::TextBox^  txtLogNumber;
	private: System::Windows::Forms::Button^  btnSearch;
	private: System::Windows::Forms::Button^  btnNew;
	private: System::Windows::Forms::Button^  btnSave;
	private: System::Windows::Forms::ComboBox^  cbLogType;
	private: System::Windows::Forms::TextBox^  txtLogDate;
	private: System::Windows::Forms::Label^  lblNumEntry;
	private: System::Windows::Forms::Label^  lblLogDate;







	private: System::Windows::Forms::Label^  lblLogType;

	private: System::Windows::Forms::RichTextBox^  rtbLogDescription;
	private: System::Windows::Forms::Label^  label4;
	private: System::Windows::Forms::TextBox^  txtLogTime;
	private: System::Windows::Forms::Label^  lblUser;
	private: System::Windows::Forms::TextBox^  txtLogUser;
	private: System::Windows::Forms::Button^  btnTransform;
	private: System::Windows::Forms::Button^  btnTransformLatex;
	private: System::Windows::Forms::Button^  btnTransformHTML;
	private: System::Windows::Forms::CheckBox^  cbxPublicInternal;
	private: System::Windows::Forms::GroupBox^  gBoxType;
	private: System::Windows::Forms::RadioButton^  rbBug;


	private: System::Windows::Forms::RadioButton^  rbLog;

	private: System::Windows::Forms::BindingNavigator^  bindingNavigator1;
	private: System::Windows::Forms::ToolStripButton^  bNavAddNewItem;
	private: System::Windows::Forms::ToolStripLabel^  bNavCountItem;


	private: System::Windows::Forms::ToolStripButton^  bindingNavigatorDeleteItem;
	private: System::Windows::Forms::ToolStripButton^  bNavMoveFirstItem;

	private: System::Windows::Forms::ToolStripButton^  bNavMovePreviousItem;

	private: System::Windows::Forms::ToolStripSeparator^  bindingNavigatorSeparator;
	private: System::Windows::Forms::ToolStripTextBox^  bNavPositionItem;

	private: System::Windows::Forms::ToolStripSeparator^  bindingNavigatorSeparator1;
	private: System::Windows::Forms::ToolStripButton^  bNavMoveNextItem;
	private: System::Windows::Forms::ToolStripButton^  bNavMoveLastItem;


	private: System::Windows::Forms::ToolStripSeparator^  bindingNavigatorSeparator2;
	private: System::Windows::Forms::RadioButton^  rBConcept;
	private: System::Windows::Forms::RadioButton^  rbToDo;
	private: System::Windows::Forms::Button^  btnCancel;
	private: System::Windows::Forms::Button^  btnRefresh;
	private: System::Windows::Forms::CheckBox^  cbxBugSolved;
	private: System::Windows::Forms::Label^  lblBugSolution;
	private: System::Windows::Forms::RichTextBox^  rtbBugSolution;


	private: System::ComponentModel::IContainer^  components;

	protected: 

	private:
		/// <summary>
		/// Erforderliche Designervariable.
		/// </summary>


#pragma region Windows Form Designer generated code
		/// <summary>
		/// Erforderliche Methode für die Designerunterstützung.
		/// Der Inhalt der Methode darf nicht mit dem Code-Editor geändert werden.
		/// </summary>
		void InitializeComponent(void)
		{
			this->components = (gcnew System::ComponentModel::Container());
			System::ComponentModel::ComponentResourceManager^  resources = (gcnew System::ComponentModel::ComponentResourceManager(Form1::typeid));
			this->txtLogNumber = (gcnew System::Windows::Forms::TextBox());
			this->btnSearch = (gcnew System::Windows::Forms::Button());
			this->btnNew = (gcnew System::Windows::Forms::Button());
			this->btnSave = (gcnew System::Windows::Forms::Button());
			this->cbLogType = (gcnew System::Windows::Forms::ComboBox());
			this->txtLogDate = (gcnew System::Windows::Forms::TextBox());
			this->lblNumEntry = (gcnew System::Windows::Forms::Label());
			this->lblLogDate = (gcnew System::Windows::Forms::Label());
			this->lblLogType = (gcnew System::Windows::Forms::Label());
			this->rtbLogDescription = (gcnew System::Windows::Forms::RichTextBox());
			this->label4 = (gcnew System::Windows::Forms::Label());
			this->txtLogTime = (gcnew System::Windows::Forms::TextBox());
			this->lblUser = (gcnew System::Windows::Forms::Label());
			this->txtLogUser = (gcnew System::Windows::Forms::TextBox());
			this->btnTransform = (gcnew System::Windows::Forms::Button());
			this->btnTransformLatex = (gcnew System::Windows::Forms::Button());
			this->btnTransformHTML = (gcnew System::Windows::Forms::Button());
			this->cbxPublicInternal = (gcnew System::Windows::Forms::CheckBox());
			this->gBoxType = (gcnew System::Windows::Forms::GroupBox());
			this->rBConcept = (gcnew System::Windows::Forms::RadioButton());
			this->rbToDo = (gcnew System::Windows::Forms::RadioButton());
			this->rbBug = (gcnew System::Windows::Forms::RadioButton());
			this->rbLog = (gcnew System::Windows::Forms::RadioButton());
			this->bindingNavigator1 = (gcnew System::Windows::Forms::BindingNavigator(this->components));
			this->bNavAddNewItem = (gcnew System::Windows::Forms::ToolStripButton());
			this->bNavCountItem = (gcnew System::Windows::Forms::ToolStripLabel());
			this->bindingNavigatorDeleteItem = (gcnew System::Windows::Forms::ToolStripButton());
			this->bNavMoveFirstItem = (gcnew System::Windows::Forms::ToolStripButton());
			this->bNavMovePreviousItem = (gcnew System::Windows::Forms::ToolStripButton());
			this->bindingNavigatorSeparator = (gcnew System::Windows::Forms::ToolStripSeparator());
			this->bNavPositionItem = (gcnew System::Windows::Forms::ToolStripTextBox());
			this->bindingNavigatorSeparator1 = (gcnew System::Windows::Forms::ToolStripSeparator());
			this->bNavMoveNextItem = (gcnew System::Windows::Forms::ToolStripButton());
			this->bNavMoveLastItem = (gcnew System::Windows::Forms::ToolStripButton());
			this->bindingNavigatorSeparator2 = (gcnew System::Windows::Forms::ToolStripSeparator());
			this->btnCancel = (gcnew System::Windows::Forms::Button());
			this->btnRefresh = (gcnew System::Windows::Forms::Button());
			this->cbxBugSolved = (gcnew System::Windows::Forms::CheckBox());
			this->lblBugSolution = (gcnew System::Windows::Forms::Label());
			this->rtbBugSolution = (gcnew System::Windows::Forms::RichTextBox());
			this->gBoxType->SuspendLayout();
			(cli::safe_cast<System::ComponentModel::ISupportInitialize^  >(this->bindingNavigator1))->BeginInit();
			this->bindingNavigator1->SuspendLayout();
			this->SuspendLayout();
			// 
			// txtLogNumber
			// 
			this->txtLogNumber->Enabled = false;
			this->txtLogNumber->Location = System::Drawing::Point(107, 101);
			this->txtLogNumber->Name = L"txtLogNumber";
			this->txtLogNumber->Size = System::Drawing::Size(216, 20);
			this->txtLogNumber->TabIndex = 0;
			// 
			// btnSearch
			// 
			this->btnSearch->Location = System::Drawing::Point(540, 61);
			this->btnSearch->Name = L"btnSearch";
			this->btnSearch->Size = System::Drawing::Size(66, 29);
			this->btnSearch->TabIndex = 1;
			this->btnSearch->Text = L"Search";
			this->btnSearch->UseVisualStyleBackColor = true;
			this->btnSearch->Visible = false;
			// 
			// btnNew
			// 
			this->btnNew->Location = System::Drawing::Point(612, 96);
			this->btnNew->Name = L"btnNew";
			this->btnNew->Size = System::Drawing::Size(66, 29);
			this->btnNew->TabIndex = 2;
			this->btnNew->Text = L"New";
			this->btnNew->UseVisualStyleBackColor = true;
			this->btnNew->Click += gcnew System::EventHandler(this, &Form1::btnNew_Click);
			// 
			// btnSave
			// 
			this->btnSave->Location = System::Drawing::Point(684, 96);
			this->btnSave->Name = L"btnSave";
			this->btnSave->Size = System::Drawing::Size(66, 29);
			this->btnSave->TabIndex = 3;
			this->btnSave->Text = L"Save";
			this->btnSave->UseVisualStyleBackColor = true;
			this->btnSave->Click += gcnew System::EventHandler(this, &Form1::btnSave_Click);
			// 
			// cbLogType
			// 
			this->cbLogType->AccessibleDescription = L"change, extension, new feature";
			this->cbLogType->AutoCompleteCustomSource->AddRange(gcnew cli::array< System::String^  >(3) {L"change", L"extension", L"new feature"});
			this->cbLogType->AutoCompleteMode = System::Windows::Forms::AutoCompleteMode::Append;
			this->cbLogType->AutoCompleteSource = System::Windows::Forms::AutoCompleteSource::ListItems;
			this->cbLogType->FormattingEnabled = true;
			this->cbLogType->Items->AddRange(gcnew cli::array< System::Object^  >(3) {L"change", L"extension", L"new feature"});
			this->cbLogType->Location = System::Drawing::Point(107, 179);
			this->cbLogType->Name = L"cbLogType";
			this->cbLogType->RightToLeft = System::Windows::Forms::RightToLeft::No;
			this->cbLogType->Size = System::Drawing::Size(216, 21);
			this->cbLogType->TabIndex = 4;
			this->cbLogType->Text = L"change";
			// 
			// txtLogDate
			// 
			this->txtLogDate->Location = System::Drawing::Point(107, 153);
			this->txtLogDate->Name = L"txtLogDate";
			this->txtLogDate->Size = System::Drawing::Size(110, 20);
			this->txtLogDate->TabIndex = 5;
			// 
			// lblNumEntry
			// 
			this->lblNumEntry->AutoSize = true;
			this->lblNumEntry->Location = System::Drawing::Point(13, 104);
			this->lblNumEntry->Name = L"lblNumEntry";
			this->lblNumEntry->Size = System::Drawing::Size(65, 13);
			this->lblNumEntry->TabIndex = 6;
			this->lblNumEntry->Text = L"Log-Number";
			// 
			// lblLogDate
			// 
			this->lblLogDate->AutoSize = true;
			this->lblLogDate->Location = System::Drawing::Point(13, 156);
			this->lblLogDate->Name = L"lblLogDate";
			this->lblLogDate->Size = System::Drawing::Size(79, 13);
			this->lblLogDate->TabIndex = 7;
			this->lblLogDate->Text = L"Log-Date/Time";
			// 
			// lblLogType
			// 
			this->lblLogType->AutoSize = true;
			this->lblLogType->Location = System::Drawing::Point(13, 182);
			this->lblLogType->Name = L"lblLogType";
			this->lblLogType->Size = System::Drawing::Size(52, 13);
			this->lblLogType->TabIndex = 8;
			this->lblLogType->Text = L"Log-Type";
			// 
			// rtbLogDescription
			// 
			this->rtbLogDescription->EnableAutoDragDrop = true;
			this->rtbLogDescription->Location = System::Drawing::Point(107, 206);
			this->rtbLogDescription->Name = L"rtbLogDescription";
			this->rtbLogDescription->ScrollBars = System::Windows::Forms::RichTextBoxScrollBars::ForcedBoth;
			this->rtbLogDescription->Size = System::Drawing::Size(643, 224);
			this->rtbLogDescription->TabIndex = 10;
			this->rtbLogDescription->Text = L"";
			// 
			// label4
			// 
			this->label4->AutoSize = true;
			this->label4->Location = System::Drawing::Point(12, 206);
			this->label4->Name = L"label4";
			this->label4->Size = System::Drawing::Size(60, 13);
			this->label4->TabIndex = 11;
			this->label4->Text = L"Description";
			// 
			// txtLogTime
			// 
			this->txtLogTime->Location = System::Drawing::Point(223, 153);
			this->txtLogTime->Name = L"txtLogTime";
			this->txtLogTime->Size = System::Drawing::Size(100, 20);
			this->txtLogTime->TabIndex = 12;
			// 
			// lblUser
			// 
			this->lblUser->AutoSize = true;
			this->lblUser->Location = System::Drawing::Point(13, 130);
			this->lblUser->Name = L"lblUser";
			this->lblUser->Size = System::Drawing::Size(29, 13);
			this->lblUser->TabIndex = 14;
			this->lblUser->Text = L"User";
			// 
			// txtLogUser
			// 
			this->txtLogUser->AcceptsTab = true;
			this->txtLogUser->Location = System::Drawing::Point(107, 127);
			this->txtLogUser->Name = L"txtLogUser";
			this->txtLogUser->Size = System::Drawing::Size(216, 20);
			this->txtLogUser->TabIndex = 13;
			// 
			// btnTransform
			// 
			this->btnTransform->Location = System::Drawing::Point(540, 131);
			this->btnTransform->Name = L"btnTransform";
			this->btnTransform->Size = System::Drawing::Size(66, 29);
			this->btnTransform->TabIndex = 15;
			this->btnTransform->Text = L"Transform";
			this->btnTransform->UseVisualStyleBackColor = true;
			this->btnTransform->Click += gcnew System::EventHandler(this, &Form1::btnTransform_Click);
			// 
			// btnTransformLatex
			// 
			this->btnTransformLatex->Location = System::Drawing::Point(612, 166);
			this->btnTransformLatex->Name = L"btnTransformLatex";
			this->btnTransformLatex->Size = System::Drawing::Size(138, 29);
			this->btnTransformLatex->TabIndex = 16;
			this->btnTransformLatex->Text = L"Transform to Latex";
			this->btnTransformLatex->UseVisualStyleBackColor = true;
			this->btnTransformLatex->Visible = false;
			this->btnTransformLatex->Click += gcnew System::EventHandler(this, &Form1::btnTransformLatex_Click);
			// 
			// btnTransformHTML
			// 
			this->btnTransformHTML->Location = System::Drawing::Point(612, 131);
			this->btnTransformHTML->Name = L"btnTransformHTML";
			this->btnTransformHTML->Size = System::Drawing::Size(138, 29);
			this->btnTransformHTML->TabIndex = 17;
			this->btnTransformHTML->Text = L"Transform to HTML";
			this->btnTransformHTML->UseVisualStyleBackColor = true;
			this->btnTransformHTML->Click += gcnew System::EventHandler(this, &Form1::btnTransformHTML_Click);
			// 
			// cbxPublicInternal
			// 
			this->cbxPublicInternal->AutoSize = true;
			this->cbxPublicInternal->Location = System::Drawing::Point(341, 181);
			this->cbxPublicInternal->Name = L"cbxPublicInternal";
			this->cbxPublicInternal->Size = System::Drawing::Size(54, 17);
			this->cbxPublicInternal->TabIndex = 18;
			this->cbxPublicInternal->Text = L"public";
			this->cbxPublicInternal->UseVisualStyleBackColor = true;
			// 
			// gBoxType
			// 
			this->gBoxType->Controls->Add(this->rBConcept);
			this->gBoxType->Controls->Add(this->rbToDo);
			this->gBoxType->Controls->Add(this->rbBug);
			this->gBoxType->Controls->Add(this->rbLog);
			this->gBoxType->Location = System::Drawing::Point(15, 27);
			this->gBoxType->Name = L"gBoxType";
			this->gBoxType->Size = System::Drawing::Size(194, 63);
			this->gBoxType->TabIndex = 19;
			this->gBoxType->TabStop = false;
			this->gBoxType->Text = L"Type";
			// 
			// rBConcept
			// 
			this->rBConcept->AutoSize = true;
			this->rBConcept->Enabled = false;
			this->rBConcept->Location = System::Drawing::Point(107, 40);
			this->rBConcept->Name = L"rBConcept";
			this->rBConcept->Size = System::Drawing::Size(65, 17);
			this->rBConcept->TabIndex = 3;
			this->rBConcept->Text = L"Concept";
			this->rBConcept->UseVisualStyleBackColor = true;
			// 
			// rbToDo
			// 
			this->rbToDo->AutoSize = true;
			this->rbToDo->Enabled = false;
			this->rbToDo->Location = System::Drawing::Point(107, 18);
			this->rbToDo->Name = L"rbToDo";
			this->rbToDo->Size = System::Drawing::Size(52, 17);
			this->rbToDo->TabIndex = 2;
			this->rbToDo->Text = L"ToDo";
			this->rbToDo->UseVisualStyleBackColor = true;
			// 
			// rbBug
			// 
			this->rbBug->AutoSize = true;
			this->rbBug->Location = System::Drawing::Point(11, 40);
			this->rbBug->Name = L"rbBug";
			this->rbBug->Size = System::Drawing::Size(44, 17);
			this->rbBug->TabIndex = 1;
			this->rbBug->Text = L"Bug";
			this->rbBug->UseVisualStyleBackColor = true;
			// 
			// rbLog
			// 
			this->rbLog->AutoSize = true;
			this->rbLog->Checked = true;
			this->rbLog->Location = System::Drawing::Point(11, 18);
			this->rbLog->Name = L"rbLog";
			this->rbLog->Size = System::Drawing::Size(43, 17);
			this->rbLog->TabIndex = 0;
			this->rbLog->TabStop = true;
			this->rbLog->Text = L"Log";
			this->rbLog->UseVisualStyleBackColor = true;
			this->rbLog->CheckedChanged += gcnew System::EventHandler(this, &Form1::rbLog_CheckedChanged);
			// 
			// bindingNavigator1
			// 
			this->bindingNavigator1->AddNewItem = this->bNavAddNewItem;
			this->bindingNavigator1->CountItem = this->bNavCountItem;
			this->bindingNavigator1->DeleteItem = this->bindingNavigatorDeleteItem;
			this->bindingNavigator1->Items->AddRange(gcnew cli::array< System::Windows::Forms::ToolStripItem^  >(11) {this->bNavMoveFirstItem, 
				this->bNavMovePreviousItem, this->bindingNavigatorSeparator, this->bNavPositionItem, this->bNavCountItem, this->bindingNavigatorSeparator1, 
				this->bNavMoveNextItem, this->bNavMoveLastItem, this->bindingNavigatorSeparator2, this->bNavAddNewItem, this->bindingNavigatorDeleteItem});
			this->bindingNavigator1->Location = System::Drawing::Point(0, 0);
			this->bindingNavigator1->MoveFirstItem = this->bNavMoveFirstItem;
			this->bindingNavigator1->MoveLastItem = this->bNavMoveLastItem;
			this->bindingNavigator1->MoveNextItem = this->bNavMoveNextItem;
			this->bindingNavigator1->MovePreviousItem = this->bNavMovePreviousItem;
			this->bindingNavigator1->Name = L"bindingNavigator1";
			this->bindingNavigator1->PositionItem = this->bNavPositionItem;
			this->bindingNavigator1->Size = System::Drawing::Size(768, 25);
			this->bindingNavigator1->TabIndex = 21;
			this->bindingNavigator1->Text = L"bindingNavigator1";
			// 
			// bNavAddNewItem
			// 
			this->bNavAddNewItem->DisplayStyle = System::Windows::Forms::ToolStripItemDisplayStyle::Image;
			this->bNavAddNewItem->Image = (cli::safe_cast<System::Drawing::Image^  >(resources->GetObject(L"bNavAddNewItem.Image")));
			this->bNavAddNewItem->Name = L"bNavAddNewItem";
			this->bNavAddNewItem->RightToLeftAutoMirrorImage = true;
			this->bNavAddNewItem->Size = System::Drawing::Size(23, 22);
			this->bNavAddNewItem->Text = L"Neu hinzufügen";
			this->bNavAddNewItem->Click += gcnew System::EventHandler(this, &Form1::bNavAddNewItem_Click);
			// 
			// bNavCountItem
			// 
			this->bNavCountItem->Name = L"bNavCountItem";
			this->bNavCountItem->Size = System::Drawing::Size(44, 22);
			this->bNavCountItem->Text = L"von {0}";
			this->bNavCountItem->ToolTipText = L"Die Gesamtanzahl der Elemente.";
			// 
			// bindingNavigatorDeleteItem
			// 
			this->bindingNavigatorDeleteItem->DisplayStyle = System::Windows::Forms::ToolStripItemDisplayStyle::Image;
			this->bindingNavigatorDeleteItem->Image = (cli::safe_cast<System::Drawing::Image^  >(resources->GetObject(L"bindingNavigatorDeleteItem.Image")));
			this->bindingNavigatorDeleteItem->Name = L"bindingNavigatorDeleteItem";
			this->bindingNavigatorDeleteItem->RightToLeftAutoMirrorImage = true;
			this->bindingNavigatorDeleteItem->Size = System::Drawing::Size(23, 22);
			this->bindingNavigatorDeleteItem->Text = L"Löschen";
			this->bindingNavigatorDeleteItem->Visible = false;
			// 
			// bNavMoveFirstItem
			// 
			this->bNavMoveFirstItem->DisplayStyle = System::Windows::Forms::ToolStripItemDisplayStyle::Image;
			this->bNavMoveFirstItem->Image = (cli::safe_cast<System::Drawing::Image^  >(resources->GetObject(L"bNavMoveFirstItem.Image")));
			this->bNavMoveFirstItem->Name = L"bNavMoveFirstItem";
			this->bNavMoveFirstItem->RightToLeftAutoMirrorImage = true;
			this->bNavMoveFirstItem->Size = System::Drawing::Size(23, 22);
			this->bNavMoveFirstItem->Text = L"Erste verschieben";
			this->bNavMoveFirstItem->Click += gcnew System::EventHandler(this, &Form1::bNavMoveFirstItem_Click);
			// 
			// bNavMovePreviousItem
			// 
			this->bNavMovePreviousItem->DisplayStyle = System::Windows::Forms::ToolStripItemDisplayStyle::Image;
			this->bNavMovePreviousItem->Image = (cli::safe_cast<System::Drawing::Image^  >(resources->GetObject(L"bNavMovePreviousItem.Image")));
			this->bNavMovePreviousItem->Name = L"bNavMovePreviousItem";
			this->bNavMovePreviousItem->RightToLeftAutoMirrorImage = true;
			this->bNavMovePreviousItem->Size = System::Drawing::Size(23, 22);
			this->bNavMovePreviousItem->Text = L"Vorherige verschieben";
			this->bNavMovePreviousItem->Click += gcnew System::EventHandler(this, &Form1::bNavMovePreviousItem_Click);
			// 
			// bindingNavigatorSeparator
			// 
			this->bindingNavigatorSeparator->Name = L"bindingNavigatorSeparator";
			this->bindingNavigatorSeparator->Size = System::Drawing::Size(6, 25);
			// 
			// bNavPositionItem
			// 
			this->bNavPositionItem->AccessibleName = L"Position";
			this->bNavPositionItem->AutoSize = false;
			this->bNavPositionItem->Name = L"bNavPositionItem";
			this->bNavPositionItem->Size = System::Drawing::Size(50, 23);
			this->bNavPositionItem->Text = L"0";
			this->bNavPositionItem->ToolTipText = L"Aktuelle Position";
			this->bNavPositionItem->TextChanged += gcnew System::EventHandler(this, &Form1::bNavPositionItem_TextChanged);
			// 
			// bindingNavigatorSeparator1
			// 
			this->bindingNavigatorSeparator1->Name = L"bindingNavigatorSeparator1";
			this->bindingNavigatorSeparator1->Size = System::Drawing::Size(6, 25);
			// 
			// bNavMoveNextItem
			// 
			this->bNavMoveNextItem->DisplayStyle = System::Windows::Forms::ToolStripItemDisplayStyle::Image;
			this->bNavMoveNextItem->Image = (cli::safe_cast<System::Drawing::Image^  >(resources->GetObject(L"bNavMoveNextItem.Image")));
			this->bNavMoveNextItem->Name = L"bNavMoveNextItem";
			this->bNavMoveNextItem->RightToLeftAutoMirrorImage = true;
			this->bNavMoveNextItem->Size = System::Drawing::Size(23, 22);
			this->bNavMoveNextItem->Text = L"Nächste verschieben";
			this->bNavMoveNextItem->Click += gcnew System::EventHandler(this, &Form1::bNavMoveNextItem_Click);
			// 
			// bNavMoveLastItem
			// 
			this->bNavMoveLastItem->DisplayStyle = System::Windows::Forms::ToolStripItemDisplayStyle::Image;
			this->bNavMoveLastItem->Image = (cli::safe_cast<System::Drawing::Image^  >(resources->GetObject(L"bNavMoveLastItem.Image")));
			this->bNavMoveLastItem->Name = L"bNavMoveLastItem";
			this->bNavMoveLastItem->RightToLeftAutoMirrorImage = true;
			this->bNavMoveLastItem->Size = System::Drawing::Size(23, 22);
			this->bNavMoveLastItem->Text = L"Letzte verschieben";
			this->bNavMoveLastItem->Click += gcnew System::EventHandler(this, &Form1::bNavMoveLastItem_Click);
			// 
			// bindingNavigatorSeparator2
			// 
			this->bindingNavigatorSeparator2->Name = L"bindingNavigatorSeparator2";
			this->bindingNavigatorSeparator2->Size = System::Drawing::Size(6, 25);
			// 
			// btnCancel
			// 
			this->btnCancel->Location = System::Drawing::Point(540, 96);
			this->btnCancel->Name = L"btnCancel";
			this->btnCancel->Size = System::Drawing::Size(66, 29);
			this->btnCancel->TabIndex = 22;
			this->btnCancel->Text = L"Cancel";
			this->btnCancel->UseVisualStyleBackColor = true;
			this->btnCancel->Click += gcnew System::EventHandler(this, &Form1::btnCancel_Click);
			// 
			// btnRefresh
			// 
			this->btnRefresh->Location = System::Drawing::Point(684, 61);
			this->btnRefresh->Name = L"btnRefresh";
			this->btnRefresh->Size = System::Drawing::Size(66, 29);
			this->btnRefresh->TabIndex = 23;
			this->btnRefresh->Text = L"Refresh";
			this->btnRefresh->UseVisualStyleBackColor = true;
			this->btnRefresh->Click += gcnew System::EventHandler(this, &Form1::btnRefresh_Click);
			// 
			// cbxBugSolved
			// 
			this->cbxBugSolved->AutoSize = true;
			this->cbxBugSolved->Location = System::Drawing::Point(341, 158);
			this->cbxBugSolved->Name = L"cbxBugSolved";
			this->cbxBugSolved->Size = System::Drawing::Size(78, 17);
			this->cbxBugSolved->TabIndex = 24;
			this->cbxBugSolved->Text = L"bug solved";
			this->cbxBugSolved->UseVisualStyleBackColor = true;
			this->cbxBugSolved->Visible = false;
			this->cbxBugSolved->CheckedChanged += gcnew System::EventHandler(this, &Form1::cbxBugSolved_CheckedChanged);
			// 
			// lblBugSolution
			// 
			this->lblBugSolution->AutoSize = true;
			this->lblBugSolution->Location = System::Drawing::Point(12, 358);
			this->lblBugSolution->Name = L"lblBugSolution";
			this->lblBugSolution->Size = System::Drawing::Size(45, 13);
			this->lblBugSolution->TabIndex = 25;
			this->lblBugSolution->Text = L"Solution";
			// 
			// rtbBugSolution
			// 
			this->rtbBugSolution->EnableAutoDragDrop = true;
			this->rtbBugSolution->Location = System::Drawing::Point(107, 358);
			this->rtbBugSolution->Name = L"rtbBugSolution";
			this->rtbBugSolution->ScrollBars = System::Windows::Forms::RichTextBoxScrollBars::ForcedBoth;
			this->rtbBugSolution->Size = System::Drawing::Size(643, 72);
			this->rtbBugSolution->TabIndex = 26;
			this->rtbBugSolution->Text = L"";
			this->rtbBugSolution->TextChanged += gcnew System::EventHandler(this, &Form1::rtbBugSolution_TextChanged);
			// 
			// Form1
			// 
			this->AutoScaleDimensions = System::Drawing::SizeF(6, 13);
			this->AutoScaleMode = System::Windows::Forms::AutoScaleMode::Font;
			this->ClientSize = System::Drawing::Size(768, 442);
			this->Controls->Add(this->rtbBugSolution);
			this->Controls->Add(this->lblBugSolution);
			this->Controls->Add(this->cbxBugSolved);
			this->Controls->Add(this->btnRefresh);
			this->Controls->Add(this->btnCancel);
			this->Controls->Add(this->bindingNavigator1);
			this->Controls->Add(this->gBoxType);
			this->Controls->Add(this->cbxPublicInternal);
			this->Controls->Add(this->btnTransformHTML);
			this->Controls->Add(this->btnTransformLatex);
			this->Controls->Add(this->btnTransform);
			this->Controls->Add(this->lblUser);
			this->Controls->Add(this->txtLogUser);
			this->Controls->Add(this->txtLogTime);
			this->Controls->Add(this->label4);
			this->Controls->Add(this->rtbLogDescription);
			this->Controls->Add(this->lblLogType);
			this->Controls->Add(this->lblLogDate);
			this->Controls->Add(this->lblNumEntry);
			this->Controls->Add(this->txtLogDate);
			this->Controls->Add(this->cbLogType);
			this->Controls->Add(this->btnSave);
			this->Controls->Add(this->btnNew);
			this->Controls->Add(this->btnSearch);
			this->Controls->Add(this->txtLogNumber);
			this->Name = L"Form1";
			this->Text = L"Form1";
			this->Load += gcnew System::EventHandler(this, &Form1::Form1_Load);
			this->gBoxType->ResumeLayout(false);
			this->gBoxType->PerformLayout();
			(cli::safe_cast<System::ComponentModel::ISupportInitialize^  >(this->bindingNavigator1))->EndInit();
			this->bindingNavigator1->ResumeLayout(false);
			this->bindingNavigator1->PerformLayout();
			this->ResumeLayout(false);
			this->PerformLayout();

		}
public:
		mystr* user;
		mystr* log_filename;
		mystr* bug_filename;
		mystr* todo_filename;
		mystr* concept_filename;
		mystr* latex_table_filename;
		mystr* html_filename;
		mystr* html_public_filename;
		mystr* html_bug_filename;
		mystr* html_bug_public_filename;
		mystr* hotint_version_filename;
		//LogData* lastLogData;
		TArrayDynamic<LogData>* logData;
		TArrayDynamic<BugData>* bugData;
		bool isLog;
		bool isBug;
		bool isNew;
		bool bugChanged;
		int curBug;
		//TArrayDynamic<mystr>* header;

private:

		void Destroy()
		{
			delete user;
			delete log_filename;
			delete bug_filename;
			delete todo_filename;
			delete concept_filename;
			delete latex_table_filename;
			delete html_filename;
			delete html_public_filename;
			delete html_bug_filename;
			delete html_bug_public_filename;
			delete hotint_version_filename;
		}


		void ExtractLog(ifstream& is, LogData& data)
		{
			mystr s1, s2, s3, s4;
			char c1, c2;
			char input[1000];
			int pos;
			//is >> str;
			//is >> line;
			//is >> line;
			is >> c1;
			is >> c2;
			//mystr tmp = str.SubString(0,1);
			//if (tmp == mystr(IDENT_NEW_LOG))
			if (c1 == '@' && c2 == '@')
			{
				//int posend = str.Find(2,',');

				//last char is a ','; the remaining part is the number
				//data.log_num =  str.SubString(2,posend-1).MakeInt();
				//is >> line;
				is >> data.ver_num;
				is >> c1;
				if (c1 == '.')
				{
					is >> data.rel_num;
					is >> c1;
					is >> data.log_num;
					is >> c1;
				}
				else
				{
					data.log_num = data.ver_num;
					data.ver_num = 1;
					data.rel_num = 1;
				}
				is >> data.date.year;
				is >> c1;
				is >> data.date.month;
				is >> c1;
				is >> data.date.day;
				//is >> s1;
				is.getline(input,1000);
				s1 = mystr(input);
				pos = 2;  //skip ", "
				s1.GetUntil(pos,',',data.time);
				++pos;
				s1.GetUntil(pos,',',data.user);
				++pos;
				s1.GetUntil(pos,',',data.type);
				data.type.EraseSpacesHeadTail();
				++pos;
				if (data.type.Length() == s1.Length())
				{
					data.entry = data.type;
					data.type = "";
					data.public_internal = false;
				}
				else
				{
					//int public_internal = 0;
					int pos_old = pos;
					s1.GetUntil(pos,',',s2);
					s2.EraseSpacesHeadTail();
					int tmpInt = s2.MakeInt();
					if (s2.Length() == 1 && (tmpInt == 0 || tmpInt == 1))
					{
						data.public_internal = tmpInt == 1 ? true : false; 
						++pos;
						//str.GetUntilEOL(pos,',',data.entry);
						data.entry = s1.SubString(pos,s1.Length()-1);
					}
					else
					{
						data.public_internal = false;
						
						//add rest of string
						if (pos < s1.Length())
						{
							data.entry = s1.SubString(pos_old, s1.Length()-1);
						}
						else
							data.entry = s2;
					}
				}
				int eoflog = is.eof();
				while(!eoflog && !is.eof())
				{
					is.getline(input,1000);
					s1 = mystr(input);
					//is >> s1;
					//s1.SetLength(0); 
					//is.getline(s1.c_str(), 256);	
					if (s1.SubString(0,1) == mystr(IDENT_NEW_LOG))
					{ eoflog = 1; }
					else
					{	data.entry += mystr("\n") + s1; 
					  //is >> line;
					}
				}
				int len = strlen(input);
				//is.seekg(-1, ios_base::cur);
				is.seekg(-len-2, ios_base::cur);
				//is.getline(input,1000);
			}
			else
			{
				data.ver_num = -1;
				data.rel_num = 0;
				data.log_num = 0;
				data.date.day = 0;
				data.date.month = 0;
				data.date.year = 0;
				data.public_internal = false;
			}
			CleanUpLog(data);
			//return -1;
		}

		int ExtractBug(ifstream& is, BugData& data, int is_solved)
		{
			int eof_open_bugs = 0;
			mystr s1, s2, s3, s4;
			char c1, c2, c3;
			char input[1000];
			int pos;
			//is >> str;
			//is >> line;
			//is >> line;
			is >> c1;
			is >> c2;
			is >> c3;
			//mystr tmp = str.SubString(0,1);
			//if (tmp == mystr(IDENT_NEW_LOG))
			if (c1 == '@' && c2 == '@' && c3 == 'B')
			{
				//int posend = str.Find(2,',');

				//last char is a ','; the remaining part is the number
				//data.log_num =  str.SubString(2,posend-1).MakeInt();
				//is >> line;
				is >> data.ver_num;
				is >> c1;
				if (c1 == '.')
				{
					is >> data.rel_num;
					is >> c1;
					is >> data.log_num;
					is >> c1;
				}
				else
				{
					data.log_num = data.ver_num;
					data.ver_num = 2;
					data.rel_num = 1;
				}
				is >> data.date.year;
				is >> c1;
				is >> data.date.month;
				is >> c1;
				is >> data.date.day;
				//is >> s1;
				is.getline(input,1000);
				s1 = mystr(input);
				pos = 2;  //skip ", "
				s1.GetUntil(pos,',',data.time);
				++pos;
				s1.GetUntil(pos,',',data.user);
				data.user.EraseSpacesHeadTail();
				++pos;
				mystr totaltmp;
				s1.GetUntil(pos,',', totaltmp);
				totaltmp.EraseSpacesHeadTail();
				if (totaltmp == mystr("0") || totaltmp == mystr("1"))
				{
					data.public_internal = (totaltmp== mystr("1")) ? true : false;
					++pos;
				}
				else
				{
					data.public_internal = false;
				}
				data.entry = s1.SubString(pos,s1.Length()-1);
				int eofbug = is.eof();
				mystr entry_addon;
				while(!eofbug && !is.eof())
				{
					is.getline(input,1000);
					s1 = mystr(input);
					//is >> s1;
					//s1.SetLength(0); 
					//is.getline(s1.c_str(), 256);	
					if (s1.SubString(0,2) == mystr(IDENT_NEW_BUG))
					{ 
						eofbug = 1; 
						int len = strlen(input);
						//is.seekg(-1, ios_base::cur);
						is.seekg(-len-2, ios_base::cur);
					}
					else if (s1 == mystr(IDENT_START_SOLVED_BUGS))
					{ 
						is.getline(input,1000);
						s1 = mystr(input);
						if ( s1 == mystr(IDENT_BUGS_ADD_SOLUTION))
						{
							eof_open_bugs = 1;
							eofbug = 1;
						}
						else
							data.entry += IDENT_START_SOLVED_BUGS + mystr("\n") + s1;
					}
					else
					{	
						//if there are additional \n add them only if another entry is added lateron
						if (s1.Length() == 0) 
						{
							entry_addon += mystr("\n");
						}
						else
							data.entry += mystr("\n") + entry_addon + mystr("\n") + s1; 
					  //is >> line;
					}
				}	
				data.solved = (bool)is_solved;
				if (is_solved)
					ExtractSolvedBug(data);
			}
			else
			{
				data.ver_num = -1;
				data.rel_num = 0;
				data.log_num = 0;
				data.date.day = 0;
				data.date.month = 0;
				data.date.year = 0;
				data.solved = false;
			}
			CleanUpBug(data);
			return eof_open_bugs;
		}

		void ExtractSolvedBug(BugData& data)
		{
			int pos = data.entry.Find(mystr(IDENT_SOLUTION));
			if (pos > 0)
			{
				data.solution = data.entry.Right(data.entry.Length() - pos - mystr(IDENT_SOLUTION).Length());
				data.solution.EraseSpacesHeadTail();
				mystr tmp = data.entry.Left(pos-1);
				data.entry = tmp;
				data.entry.EraseSpacesHeadTail();
			}
		}

		mystr BuildDate(const Date& date)
		{
			mystr mydate(date.year);
			mydate += "-";
			if (date.month < 10)
				mydate += "0";
			mydate += date.month;
			mydate += "-";
			if (date.day < 10)
				mydate += "0";
			mydate += date.day;
			return mydate;
		}

		void GetDate(mystr& date, const Date& mydate)
		{
			date = BuildDate(mydate);
		}

		//used for both logs and bugs (log_num might be the acutal log-number or the bug-number)
		mystr BuildLogNum(const BaseData& data)
		{
			mystr rv;
			rv = mystr(data.ver_num) + "." + mystr(data.rel_num) + "." + mystr(data.log_num);
			return rv;
		}
		
		/*mystr BuildBugNum(const BugData& data)
		{
			mystr rv;
			rv = mystr(data.ver_num) + "." + mystr(data.rel_num) + "." + mystr(data.bug_num);
			return rv;
		}*/

		void FillData(LogData& data)
		{
			mystr date = BuildDate(data.date);
			//System::String^ str;
			//str = date.c_str();
			//std::string ansi = StringConvA(s);
			//this->txtLogDate->Text = gcnew String(reinterpret_cast<const char*>(date.c_str()));;
			this->txtLogNumber->Text = gcnew String(BuildLogNum(data).c_str());
			this->txtLogUser->Text = gcnew String(data.user.c_str());
			this->txtLogDate->Text = gcnew String(date.c_str());
			this->txtLogTime->Text = gcnew String(data.time.c_str());
			this->cbLogType->Text = gcnew String(data.type.c_str());
			this->rtbLogDescription->Text = gcnew String(data.entry.c_str());
			this->cbxPublicInternal->Checked = data.public_internal;
			EnableFields(false);
			EnableButtons(true);
			this->btnSave->Enabled = false;
		}

		void EnableFields(bool val)
		{
			//this->txtLogNumber->Enabled = val;
			this->txtLogUser->Enabled = val;
			this->txtLogDate->Enabled = val;
			this->txtLogTime->Enabled = val;
			this->cbLogType->Enabled = val;
			//this->cbLogType->Lock
			this->rtbLogDescription->ReadOnly = !val;
			this->cbxPublicInternal->Enabled = val;
			this->bNavPositionItem->Enabled = !val;
		}

		void EnableButtons(bool val)
		{
			//this->txtLogNumber->Enabled = val;
			this->btnNew->Enabled = val;
			this->btnRefresh->Enabled = val;
			this->bNavAddNewItem->Enabled = val;
			this->btnSave->Enabled = false;
			this->btnSearch->Enabled = val;
			//this->btnTransformHTML->Enabled = false;
			//this->btnTransformHTML->Enabled = true;
			this->btnTransformHTML->Enabled = val;
			this->btnTransformLatex->Enabled = val;
			this->btnTransform->Enabled = !val;
			this->btnCancel->Enabled = false;
		}

		void EnableBNav(bool val)
		{
			this->bNavMoveFirstItem->Enabled = val;
			this->bNavMoveLastItem->Enabled = val;
			this->bNavMovePreviousItem->Enabled = val;
			this->bNavMoveNextItem->Enabled = val;
			this->bNavPositionItem->Enabled = val;
		}

		void EmptyFields()
		{
			this->txtLogNumber->Text = "";
			this->txtLogUser->Text = "";
			this->txtLogDate->Text = "";
			this->txtLogTime->Text = "";
			this->cbLogType->Text = "";
			this->rtbLogDescription->Text = "";
			this->rtbBugSolution->Text = "";
			//this->cbxInternalPublic->CheckState = Windows::Forms::CheckState::Indeterminate;
			this->cbxPublicInternal->Checked = false;
		}

		int LoadUser()
		{
			int rv = 0;
			ifstream is_user(mystr(USERNAME_FILENAME).c_str());
			//mystr str_user;
			char tmp[1000];
			if (is_user.is_open())
			{
				//is_user >> str_user;
				is_user.getline(tmp,1000);
				user = new mystr(tmp);
				//is_user >> str_user;
				is_user.getline(tmp,1000);
				latex_table_filename = new mystr(tmp);
				is_user.getline(tmp,1000);
				html_filename = new mystr(tmp);
				is_user.getline(tmp,1000);
				html_public_filename = new mystr(tmp);
				is_user.getline(tmp,1000);
				html_bug_filename = new mystr(tmp);
				is_user.getline(tmp,1000);
				html_bug_public_filename = new mystr(tmp);
				is_user.close();
			}
			else
			{
				rv = 1;
				MessageBox::Show("The user-file could not be opened!","Problem",MessageBoxButtons::OK,MessageBoxIcon::Error);
			}
			return rv;
		}

		int LoadLogFilenames()
		{
			int rv = 0;
			ifstream is_log(mystr(LOG_FILENAME).c_str());
			//mystr str_user;
			char tmp[1000];
			if (is_log.is_open())
			{
				//is_user >> str_user;
				is_log.getline(tmp,1000);
				log_filename = new mystr(tmp);
				if (is_log.eof())
				{
					rv = 1;
					MessageBox::Show("The log-config-file ended before all filenames could be determined!","Problem",MessageBoxButtons::OK,MessageBoxIcon::Error);
					is_log.close();
					return rv;
				}
				is_log.getline(tmp,1000);
				bug_filename = new mystr(tmp);
				if (is_log.eof())
				{
					rv = 1;
					MessageBox::Show("The log-config-file ended before all filenames could be determined!","Problem",MessageBoxButtons::OK,MessageBoxIcon::Error);
					is_log.close();
					return rv;
				}
				is_log.getline(tmp,1000);
				todo_filename = new mystr(tmp);
				if (is_log.eof())
				{
					rv = 1;
					MessageBox::Show("The log-config-file ended before all filenames could be determined!","Problem",MessageBoxButtons::OK,MessageBoxIcon::Error);
					is_log.close();
					return rv;
				}
				is_log.getline(tmp,1000);
				concept_filename = new mystr(tmp);
				if (is_log.eof())
				{
					rv = 1;
					MessageBox::Show("The log-config-file ended before all filenames could be determined!","Problem",MessageBoxButtons::OK,MessageBoxIcon::Error);
					is_log.close();
					return rv;
				}
				is_log.getline(tmp,1000);
				hotint_version_filename = new mystr(tmp);

				//is_user >> str_user;
				is_log.close();
			}
			else
			{
				rv = 1;
				MessageBox::Show("The log-config-file could not be opened!","Problem",MessageBoxButtons::OK,MessageBoxIcon::Error);
			}
			return rv;
		}

		int LoadLastLog()
		{
			int rv = 0;
			ShowLog(logData->Last());
			return rv;
		}

		int LoadLastBug()
		{
			int rv = 0;
			ShowBug(bugData->Last());
			return rv;
		}


		void ShowLog(LogData& logData)
		{
			CleanUpLog(logData);
			FillData(logData);
			if (logData.log_num > 0)
				EnableButtons(true);
			else
				EnableButtons(false);
			isNew = false;
		}

		void ShowBug(BugData& bugData)
		{
			CleanUpBug(bugData);
			FillBugData(bugData);
			if (bugData.log_num > 0)
			{
				EnableButtons(true);
				//bugChanged = false;
				BugChanged(false);
			}
			else
				EnableButtons(false);
			isNew = false;
		}

		void FillBugData(BugData& data)
		{
			mystr date = BuildDate(data.date);
			//System::String^ str;
			//str = date.c_str();
			//std::string ansi = StringConvA(s);
			//this->txtLogDate->Text = gcnew String(reinterpret_cast<const char*>(date.c_str()));;
			this->txtLogNumber->Text = gcnew String(BuildLogNum(data).c_str());
			this->txtLogUser->Text = gcnew String(data.user.c_str());
			this->txtLogDate->Text = gcnew String(date.c_str());
			this->txtLogTime->Text = gcnew String(data.time.c_str());
			//this->cbLogType->Text = gcnew String(data.type.c_str());
			this->rtbLogDescription->Text = gcnew String(data.entry.c_str());
			this->rtbBugSolution->Text = gcnew String(data.solution.c_str());
			this->cbxPublicInternal->Checked = data.public_internal;
			this->cbxBugSolved->Checked = data.solved;
			EnableFields(false);
			EnableButtons(true);
			this->btnSave->Enabled = false;
			
			this->cbxBugSolved->Enabled = !data.solved;
			this->rtbBugSolution->ReadOnly = data.solved;
		}


		//identType determines the type to activate (logs (1), bugs (2))
		int Initialize(int identType)
		{
			int rv = LoadUser();
			rv += LoadLogFilenames();
			if (rv)
				return rv;
			//rv += InitLogs();
			isLog = (identType == 1);
			isBug = (identType == 2);
			isNew = false;
			if (identType == 1)
			{	rv += ActivateLogs(); }
			else if (identType == 2)
			{ rv += ActivateBugs(); }
			return rv;
		}

		int InitializeLogDataNavigation()
		{
			int rv = 0;
			int num_entries = (isLog) ? logData->Length() : ((isBug) ? bugData->Length() : 0 );
			this->bNavCountItem->Text = "von {" + System::Convert::ToString(num_entries) + "}";;
			if (isLog)
			{
				ShowLog(logData->Last());
				this->bNavPositionItem->Text = System::Convert::ToString(logData->Length());
			}
			else if (isBug)
			{
				//ShowBug(bugData->Last());
				curBug = bugData->Length();
				this->bNavPositionItem->Text = System::Convert::ToString(bugData->Length());
			}
			//this->btnSave->Enabled = false;
			EnableBNav(true);
			this->bNavMoveLastItem->Enabled = false; 
			return rv;
		}

		int InitLogs()
		{
			int rv = 0;
			int num_logs = DetermineNumOfLogs();
			if (num_logs > 0)
			{
				if (logData)
					delete logData;
				logData = new TArrayDynamic<LogData>(num_logs);
				rv += ReadLogs(logData);
				rv += LoadLastLog();
				rv += InitializeLogDataNavigation();
			}
			else
			{
				EnableBNav(false);
				EnableButtons(false);
				EnableFields(false);
				this->btnTransform->Enabled = true;
			}
			return rv;
		}

		void ConvertStringText(mystr& dest, const System::String^ text)
		{
			pin_ptr<const wchar_t> str1 = PtrToStringChars(text);
			int len = wcslen(str1);
			//dest.SetLength(len);
			dest.ReSize(len);
			wcstombs(dest.c_str(), str1, len);
			dest.c_str()[len] = '\0';
		}

		void GetDate(Date& date, const mystr& str)
		{
			int pos = 0;
			mystr num;
			str.GetUntil(pos,'-',num);
			date.year = num.MakeInt();
			str.GetUntil(++pos,'-',num);
			date.month = num.MakeInt();
			num = str.SubString(++pos,str.Length()-1);
			date.day = num.MakeInt();
		}

		mystr GetActualDate()
		{
			mystr date;
	    time_t t = time(0);   // get time now
			struct tm * now = localtime( & t );
			date = mystr((now->tm_year + 1900));
			date += '-';
			if (now->tm_mon < 9)
				date += '0';
			date += mystr(now->tm_mon+1) + mystr('-');
			if (now->tm_mday < 9)
				date += '0';
			date += mystr(now->tm_mday);
			return date;
		}

		mystr GetActualTime()
		{
			mystr strtime;
	    time_t t = time(0);   // get time now
			struct tm * now = localtime( & t );
  		strtime = "";
			if (now->tm_hour < 10)
				strtime += '0';
			strtime += mystr((now->tm_hour));
			strtime += ':';
			if (now->tm_min < 10)
				strtime += '0';
			strtime += mystr(now->tm_min) + mystr(':');
			if (now->tm_sec < 10)
				strtime += '0';
			strtime += mystr(now->tm_sec);
			return strtime;
		}

		void BuildEntry(mystr& entry, const LogData& data)
		{

			entry = IDENT_NEW_LOG;
			mystr date;
			mystr fill(", ");
			GetDate(date, data.date);
			entry += BuildLogNum(data) + fill + date + fill;
			entry += data.time + fill + data.user;
			entry += fill + data.type + fill;
			entry += data.public_internal ? mystr(1) : mystr(0);
			entry += fill + data.entry;
		}

		void BuildEntry(mystr& entry, const BugData& data)
		{
			entry = IDENT_NEW_BUG;
			mystr date;
			mystr fill(", ");
			GetDate(date, data.date);
			entry += BuildLogNum(data) + fill + date + fill;
			entry += data.time + fill + data.user + fill;
			entry += data.public_internal ? mystr(1) : mystr(0);
			entry += fill + data.entry;
			if (data.solved)
			{
				entry += mystr(" ") + mystr(IDENT_SOLUTION) + mystr(" ") + data.solution;
			}
		}

		int SaveEntry(BaseData& data)
		{
			//TArrayDynamic<mystr> header2(10);
			int rv = 0;
			//int opt = ios::out | ios::app;
			//int opt = ios::in | ios::out;
			//int opt = ios::out;
			//ofstream os(mystr(CHANGES_LOG_FILENAME).c_str(), opt);
			//ifstream is(mystr(CHANGES_LOG_FILENAME).c_str(), opt);
			//ofstream os(log_filename->c_str(), opt);
			//ifstream is(log_filename->c_str(), opt);
			//mystr entry;
			//
			//BuildEntry(entry, data);
			////char tmp[1000];
			////mystr line;
			//if (os.is_open() && is.is_open())
			//{

			//	////os.tellp()
			//	//os.seekp(100,ios_base::beg);
			//	//os.put('Z');

			//	////os << entry;
			//	//os.close();
			//	//return 1;
			//	int stop = 0;
			//	int i = 0;
			//	while (!stop)
			//	{
			//		++i;
			//		is.getline(tmp, 1000);
			//		line = tmp;
			//		//os.rdbuf()->
			//		
			//		//if (line.SubString(0,1) == mystr(IDENT_NEW_LOG))
			//		if (line == mystr(IDENT_START_LOG))
			//		{
			//			stop++;
			//		}
			//		else
			//			os << line << "\n";
			//	}
			//	os << line << "\n";
			//	os << entry << "\n"; 
			//	while (!is.eof())
			//	{
			//		is.getline(tmp, 1000);
			//		line = tmp;
			//		os << line << "\n";
			//	}
			//	//filebuf* buf = os.rdbuf();
			//	//buf->sputn(entry.c_str(), entry.Length());
			//	//os.flush();
			//	//os.rdbuf()->seekpos(10);
			//	//os << entry;
			//	//os.write(entry.c_str(), entry.Length());
			//	is.close();
			//	os.close();
			//}
			//else
			//{
			//	MessageBox::Show("The log-file could not be opened!","Problem",MessageBoxButtons::OK,MessageBoxIcon::Error);
			//	rv = 1;
			//}
			//if (!rv)
			TArrayDynamic<mystr> header(5);
			if (isLog)
			{
				//LogData& tmpLog = dynamic_cast<LogData&> (data);
				logData->Add((LogData&)data);
				rv += ReadHeader(header);
				rv += WriteLogs(logData, header);
				rv += LoadLastLog();
			}
			else if (isBug)
			{
				//const BugData& tmpBug = dynamic_cast<const BugData&> (data);
				bugData->Add((BugData&)data);
				rv += ReadBugHeader(header);
				rv += WriteBugs(bugData, header);
				rv += LoadLastBug();
			}
			return rv;
		}

		int SaveChangedBugEntry(BugData& data)
		{
			int rv = 0;
			TArrayDynamic<mystr> header(5);
			//change the according entities (solution and solved-marker)
			(*bugData)(curBug).solution = data.solution;
			(*bugData)(curBug).solved = data.solved;
			rv += ReadBugHeader(header);
			rv += WriteBugs(bugData, header);

			ShowBug((*bugData)(curBug));
			//rv += LoadLastBug();
			return rv;
		}


		int TransformStep1()
		{
			TArrayDynamic<mystr>* data;
			int rv = 0;
			int opt = ios::in | ios::out;
			
			ifstream is(log_filename->c_str(), opt);
			char tmp[1000];
			mystr line;
			if (is.is_open())
			{
				int nlines = 0;
				while(!is.eof())
				{
					is.getline(tmp, 1000);
					nlines++;
				}
				is.seekg(0, ios_base::beg);
				//is.seekg(0, ios_base::beg);
				is.seekg(0, ios_base::beg);
				data = new TArrayDynamic<mystr>(nlines);
				while(!is.eof())
				{
					is.getline(tmp, 1000);
					line = tmp;
					data->Add(line);
				}
				is.close();
			}
			else
			{
				MessageBox::Show("Problem in transformation: The log-file could not be opened!","Problem",MessageBoxButtons::OK,MessageBoxIcon::Error);
				rv = 1;
				return rv;
			}
			ofstream os(log_filename->c_str(), opt);
			if (os.is_open()) // && is.is_open())
			{
				int stop = 0;
				int i = 0;
				int pos, act_num;
				mystr num;
				while (!stop)
				{
					++i;
					//is.getline(tmp, 1000);
					//line = tmp;
					line = (*data)(i);
					//os.rdbuf()->
					if (line == mystr(IDENT_START_LOG))
					{
						stop++;
					}
					os << line << "\n";
				}
				++i;
				int act_log = 0;
				for (; i <= data->Length(); ++i)
				{
					line = (*data)(i);
					pos = 0;
					line.GetUntil(pos,',',num);
					act_num = num.MakeInt();
					if (act_num <= 0 || ((act_log-1 != act_num) && (act_log != 0)))
						os << line << "\n";
					else
					{
						os << IDENT_NEW_LOG << "1.1." << line << "\n";
						act_log = act_num;	
					}
				}
				os.close();
			}
			else
			{
				MessageBox::Show("Problem in transformation: The log-file could not be opened!","Problem",MessageBoxButtons::OK,MessageBoxIcon::Error);
				rv = 1;
			}
			return rv;
		}

		int ReadHeader(TArrayDynamic<mystr>& header)
		{
			int rv = 0;
			int opt = ios::in;	
			ifstream is(log_filename->c_str(), opt);
			char tmp[1000];
			mystr line;
			if (is.is_open())
			{
				int stop = 0;
				int i = 0;
				while (!is.eof() && !stop && i < header.Size())
				{
					++i;
					//is.getline(tmp, 1000);
					//line = tmp;
					is.getline(tmp,1000);
					line = tmp;
					header.Add(tmp);
					//os.rdbuf()->
					if (line == mystr(IDENT_START_LOG))
						stop++;
				}
				if (stop == 0)
				{
					rv = 1;
					MessageBox::Show("Problem in reading header: START identification could not be found!","Problem",MessageBoxButtons::OK,MessageBoxIcon::Error);
				}
				is.close();
			}
			else
			{
				MessageBox::Show("Problem in reading header: The log-file could not be opened!","Problem",MessageBoxButtons::OK,MessageBoxIcon::Error);
				rv = 1;
			}
			return rv;
		}

		int ReadBugHeader(TArrayDynamic<mystr>& header)
		{
			int rv = 0;
			int opt = ios::in;	
			ifstream is(bug_filename->c_str(), opt);
			char tmp[1000];
			mystr line;
			if (is.is_open())
			{
				int stop = 0;
				int i = 0;
				while (!is.eof() && !stop && i < header.Size())
				{
					++i;
					//is.getline(tmp, 1000);
					//line = tmp;
					is.getline(tmp,1000);
					line = tmp;
					header.Add(tmp);
					//os.rdbuf()->
					if (line == mystr(IDENT_START_BUG))
						stop++;
				}
				if (stop == 0)
				{
					rv = 1;
					MessageBox::Show("Problem in reading header: START identification could not be found!","Problem",MessageBoxButtons::OK,MessageBoxIcon::Error);
				}

				is.close();
			}
			else
			{
				MessageBox::Show("Problem in reading header: The bug-file could not be opened!","Problem",MessageBoxButtons::OK,MessageBoxIcon::Error);
				rv = 1;
			}
			return rv;
		}

		int WriteLogs(TArrayDynamic<LogData>* logdata, TArrayDynamic<mystr>& header)
		{
			int rv = 0;
			int opt = ios::out;	
			ofstream os(log_filename->c_str(), opt);
			//char tmp[1000];
			mystr line, entry;
			if (os.is_open())
			{
				int stop = 0;
				int i = 0;
				for (int i = 1; i <= header.Length(); ++i)
					os << header(i) << endl;
				for (int i = logdata->Length(); i >= 1; --i)
				{
					BuildEntry(entry, (*logdata)(i));
					os << entry << endl;
				}
				os.close();
			}
			else
			{
				MessageBox::Show("Problem in transforming log-file: The log-file could not be opened!","Problem",MessageBoxButtons::OK,MessageBoxIcon::Error);
				rv = 1;
			}
			return rv;
		}

		int WriteBugs(TArrayDynamic<BugData>* bugdata, TArrayDynamic<mystr>& header)
		{
			int rv = 0;
			ofstream os(bug_filename->c_str(), ios::out);
			//char tmp[1000];
			mystr line, entry;
			if (os.is_open())
			{
				int stop = 0;
				int i = 0;
				for (int i = 1; i <= header.Length(); ++i)
					os << header(i) << endl;
				for (int i = bugdata->Length(); i >= 1; --i)
				{
					if (!(*bugdata)(i).solved)
					{
						BuildEntry(entry, (*bugdata)(i));
						os << entry << endl;
					}
				}
				os << endl << endl << endl;
				os << IDENT_START_SOLVED_BUGS << endl << IDENT_BUGS_ADD_SOLUTION << endl; 

				for (int i = bugdata->Length(); i >= 1; --i)
				{
					if ((*bugdata)(i).solved)
					{
						BuildEntry(entry, (*bugdata)(i));
						os << entry << endl;
					}
				}
				os.close();
			}
			else
			{
				MessageBox::Show("Problem in transforming bug-file: The bug-file could not be opened!","Problem",MessageBoxButtons::OK,MessageBoxIcon::Error);
				rv = 1;
			}
			return rv;
		}



		int TransformStep2()
		{
			int rv = 0;
			int num_logs = DetermineNumOfLogs();
			logData = new TArrayDynamic<LogData>(num_logs);
			rv = ReadLogs(logData);
			TArrayDynamic<mystr> header(5);
			rv += ReadHeader(header);
			rv += WriteLogs(logData, header);
			rv += LoadLastLog();
			return rv;
		}

		void CleanUpLog(LogData& data)
		{
			data.entry.EraseSpacesHeadTail();
			data.time.EraseSpacesHeadTail();
			data.user.EraseSpacesHeadTail();
			data.type.EraseSpacesHeadTail();
		}

		void CleanUpBug(BugData& data)
		{
			data.entry.EraseSpacesHeadTail();
			data.time.EraseSpacesHeadTail();
			data.user.EraseSpacesHeadTail();
			data.solution.EraseSpacesHeadTail();
		}

		void ExtractLogNum(BaseData& data, const mystr& dest)
		{
			mystr tmp;
			int pos = 0;
			dest.GetUntil(pos,'.',tmp);
			data.ver_num = tmp.MakeInt();
			++pos;
			dest.GetUntil(pos,'.',tmp);
			data.rel_num = tmp.MakeInt();
			++pos;
			tmp = dest.SubString(pos,dest.Length()-1);
			data.log_num = tmp.MakeInt();
		}


#pragma endregion
		private: 
			System::Void Form1_Load(System::Object^  sender, System::EventArgs^  e) 
			{
				Initialize(1);
			}
private: 
	System::Void btnNew_Click(System::Object^  sender, System::EventArgs^  e) 
	{
		EnableFields(true);
		EmptyFields();
		

		this->txtLogUser->Text = gcnew String(user->c_str());
		this->txtLogDate->Text = gcnew String(GetActualDate().c_str());
		this->txtLogTime->Text = gcnew String(GetActualTime().c_str());
		
		//this->txtLogNumber->Enabled = false;
		isNew = true;
		this->btnNew->Enabled = false;
		this->bNavAddNewItem->Enabled = false;
		this->btnTransformHTML->Enabled = false;
		this->btnTransformLatex->Enabled = false;
		this->btnSearch->Enabled = false;
		this->btnRefresh->Enabled = false;
		this->btnSave->Enabled = true;
		this->btnCancel->Enabled = true;
		EnableBNav(false);
		if (isLog)
		{
			LogData log(logData->Last());
			log.log_num++;
			//this->txtLogNumber->Text = gcnew String(mystr(lastLogData->log_num+1).c_str());
			this->txtLogNumber->Text = gcnew String(BuildLogNum(log).c_str());
			this->cbLogType->Text = gcnew String("change");
			this->bNavPositionItem->Text = System::Convert::ToString(logData->Length()+1);
		} 
		else if (isBug)
		{
			BugData bug(bugData->Last());
			bug.log_num++;
			//this->txtLogNumber->Text = gcnew String(mystr(lastLogData->log_num+1).c_str());
			this->txtLogNumber->Text = gcnew String(BuildLogNum(bug).c_str());
			this->bNavPositionItem->Text = System::Convert::ToString(bugData->Length()+1);
			this->cbxBugSolved->Checked = false;
			this->cbxBugSolved->Enabled = true;
			this->rtbBugSolution->ReadOnly = true;
		}
	}

	System::Void btnSave_Click(System::Object^  sender, System::EventArgs^  e) 
	{
		if (!isNew && !bugChanged)
		{
			return;
		}
		if (isLog)
		{
			SaveNewLog();
		}
		else if (isBug)
		{
			SaveBug();
		}
	}

	void SaveNewLog()
	{
		LogData data;
		mystr dest;
		ConvertStringText(dest, this->txtLogNumber->Text);
		//data.log_num = mystr(dest).MakeInt();
		ExtractLogNum(data, dest);
		ConvertStringText(dest, this->txtLogDate->Text);
		GetDate(data.date, dest);
		ConvertStringText(data.time, this->txtLogTime->Text);
		ConvertStringText(data.user, this->txtLogUser->Text);
		ConvertStringText(data.type, this->cbLogType->Text);
		data.public_internal = this->cbxPublicInternal->Checked;
		GetLogOrBugDescription(data.entry);
		CleanUpLog(data);
		int rv = SaveEntry(data);
		if (rv == 0)
		{
			isNew = false;
			EnableFields(false);
			EnableBNav(true);
			this->bNavPositionItem->Text = System::Convert::ToString(logData->Length());
		}
		rv += LoadLastLog();
		//if everything went fine enable buttons again
		if (!rv)
		{
		/* JUST TEST ... SEEMS TO BE REDUNDANT
		,this->btnNew->Enabled = true;
			this->bNavAddNewItem->Enabled = true;
			this->btnRefresh->Enabled = true;
			this->btnTransformHTML->Enabled = false;
			this->btnSave->Enabled = false;
			this->btnCancel->Enabled = false;
			this->btnTransformLatex->Enabled  = true;
			this->btnSearch->Enabled = true;*/
			rv += InitializeLogDataNavigation();
		}
		if (!rv)
		{
			rv += UpdateHotintVersion();
			MessageBox::Show("Log saved successfully!","Save",MessageBoxButtons::OK, MessageBoxIcon::Information);
		}
		else
			MessageBox::Show("Problem while saving logs!", "Save", MessageBoxButtons::OK, MessageBoxIcon::Error);
	

		rv += TransformToHTML();
	}

	
	void SaveBug()
	{
		BugData data;
		mystr dest;
		if (isNew)
		{
			ConvertStringText(dest, this->txtLogNumber->Text);
			//data.log_num = mystr(dest).MakeInt();
			ExtractLogNum(data, dest);
			ConvertStringText(dest, this->txtLogDate->Text);
			GetDate(data.date, dest);
			ConvertStringText(data.time, this->txtLogTime->Text);
			ConvertStringText(data.user, this->txtLogUser->Text);
			//ConvertStringText(data.type, this->cbLogType->Text);
			data.public_internal = this->cbxPublicInternal->Checked;
			GetLogOrBugDescription(data.entry);
		}
		data.solved = this->cbxBugSolved->Checked;
		if (data.solved)
			GetBugSolution(data.solution);

		CleanUpBug(data);
		int rv = 0;
		if (isNew)
		{
			rv = SaveEntry(data);
			if (rv == 0)
			{
				isNew = false;
				EnableFields(false);
				EnableBNav(true);
				this->bNavPositionItem->Text = System::Convert::ToString(logData->Length());
				rv += InitializeLogDataNavigation();

			}
		}
		else //bug changed
		{
			rv = SaveChangedBugEntry(data);
		}
		

		//if everything went fine enable buttons again

		if (!rv)
		{
			MessageBox::Show("Bug saved successfully!","Save",MessageBoxButtons::OK, MessageBoxIcon::Information);
		}
		else
			MessageBox::Show("Problem while saving bugs!", "Save", MessageBoxButtons::OK, MessageBoxIcon::Error);
	

		rv += TransformToHTML();
	}


	void GetLogOrBugDescription(mystr& entry)
	{
		mystr dest;

		for(int i = 0; i < this->rtbLogDescription->Lines->Length; ++i)
		{
			ConvertStringText(dest, this->rtbLogDescription->Lines[i]);
			if (i == 0)
			{	entry = dest; }
			else
			{
				entry += "\n";
				entry += dest;
			}
		}
	}

	void GetBugSolution(mystr& entry)
	{
		mystr dest;

		for(int i = 0; i < this->rtbBugSolution->Lines->Length; ++i)
		{
			ConvertStringText(dest, this->rtbBugSolution->Lines[i]);
			if (i == 0)
			{	entry = dest; }
			else
			{
				entry += "\n";
				entry += dest;
			}
		}
	}

	int UpdateHotintVersion()
	{
		int rv = 0;
		int opt = ios::out;	
		ofstream os(hotint_version_filename->c_str(), opt);
		mystr line;
		os << "//this file is autogenerated by the log-tool, don't change it manually.\n\n";
		os << "const HotintVersionInfo hotint_version(" << logData->Last().ver_num << "," << logData->Last().rel_num << "," << logData->Last().log_num << ");";
		os.close();
		return rv;
	}

	int DetermineNumOfLogs()
	{
		//determine the number of logs -> otherwise problems might occur with memory alocation
		int opt = ios::in;	
		ifstream is(log_filename->c_str(), opt);
		char tmp[1000];
		mystr line;
		int num_logs = 0;
		if (is.is_open())
		{
			int stop = 0;
			while (!stop)
			{
				is.getline(tmp,1000);
				line = tmp;
				if (line == mystr(IDENT_START_LOG))
					stop++;
			}
			while(!is.eof())
			{
				is.getline(tmp, 1000);
				line = tmp;
				if (line.SubString(0,1) == mystr(IDENT_NEW_LOG))
					num_logs++;
			}
		}
		else
		{
				MessageBox::Show("Could not open log-filename. Probably, the log-path is wrong!","DetermineNumOfLogs",MessageBoxButtons::OK,MessageBoxIcon::Error);

		}
		return num_logs;
	}

	int DetermineNumOfBugs()
	{
		//determine the number of logs -> otherwise problems might occur with memory alocation
		int opt = ios::in;	
		ifstream is(bug_filename->c_str(), opt);
		char tmp[1000];
		mystr line;
		int num_bugs = 0;
		if (is.is_open())
		{
			int stop = 0;
			while (!stop)
			{
				is.getline(tmp,1000);
				line = tmp;
				if (line == mystr(IDENT_START_BUG) || is.eof())
					stop++;
			}
			while(!is.eof())
			{
				is.getline(tmp, 1000);
				line = tmp;
				if (line.SubString(0,2) == mystr(IDENT_NEW_BUG))
					num_bugs++;
			}
		}
		else
		{
				MessageBox::Show("Could not open bug-filename. Probably, the bug-path is wrong!","DetermineNumOfBugs",MessageBoxButtons::OK,MessageBoxIcon::Error);

		}
		return num_bugs;
	}


	int ReadLogs(TArrayDynamic<LogData>* logdata)
	{
		int rv = 0;
		int opt = ios::in;	
		ifstream is(log_filename->c_str(), opt);
		char tmp[1000];
		mystr line;
		if (is.is_open())
		{
			int num_logs = 0;
			if (logdata == 0)
			{
				if (num_logs == 0)
				{
					num_logs = DetermineNumOfLogs();
					logdata = new TArrayDynamic<LogData>(num_logs);
				}
			}
			else
				num_logs = logdata->Size();
			//go to start position
			int stop = 0;
			while (!stop)
			{
				is.getline(tmp,1000);
				line = tmp;
				if (line == mystr(IDENT_START_LOG))
					stop++;
			}
			//read logs
			for (int i = num_logs; i >= 1; --i)
			{
				ExtractLog(is, (*logdata)(i));
			}
			is.close();
		}
		else
		{
			MessageBox::Show("Problem in reading the log data: The log-file could not be opened!","Problem",MessageBoxButtons::OK,MessageBoxIcon::Error);
			rv = 1;
			return rv;
		}
		return rv;
	}

	int ReadBugs(TArrayDynamic<BugData>* bugdata)
	{
		int rv = 0;
		int opt = ios::in;	
		ifstream is(bug_filename->c_str(), opt);
		char tmp[1000];
		mystr line;
		if (is.is_open())
		{
			int num_bugs = 0;
			if (bugdata == 0)
			{
				if (num_bugs == 0)
				{
					num_bugs = DetermineNumOfBugs();
					bugdata = new TArrayDynamic<BugData>(num_bugs);
				}
			}
			else
				num_bugs = bugdata->Size();
			//go to start position
			int stop = 0;
			while (!stop)
			{
				is.getline(tmp,1000);
				line = tmp;
				if (line == mystr(IDENT_START_BUG))
					stop++;
			}
			//read bugs
			int eof_open_bugs; //determines if the bugs are not solved so far
			int is_solved = 0;
			for (int i = num_bugs; i >= 1; --i)
			{
				eof_open_bugs = ExtractBug(is, (*bugdata)(i), is_solved);
			  is_solved += eof_open_bugs;
			}
			is.close();
		}
		else
		{
			MessageBox::Show("Problem in reading the bug data: The bug-file could not be opened!","Problem",MessageBoxButtons::OK,MessageBoxIcon::Error);
			rv = 1;
			return rv;
		}
		rv += SortBugs();
		return rv;
	}

	//sort the bugs in ascending order of bug-numbers (since they are mixed in the unsolved and solved part
	int SortBugs()
	{
		int rv = 0;
		TArray<int> bugnums(bugData->Length());		
		TArray<int> bugindex(bugData->Length());
		for (int i=1; i<=bugData->Length(); ++i)
		{
			bugnums.Add((*bugData)(i).log_num);
			bugindex.Add(i);
		}
		Quicksort(bugnums, *bugData);
		/*TArrayDynamic<BugData> tmpBugData(*bugData);
		for (int i = 1; i <= bugindex.Length(); ++i)
		{ (*bugData)(i) = tmpBugData(bugindex(i)); }		*/
		return rv;
	}

	mystr AdaptStringToLatex(const mystr& str)
	{
		mystr rv(str);
		//rv.Replace("\","\\ textbackslash ");
		rv.Replace("\\","\\textbackslash " );
		rv.Replace("&", "\\&");
		rv.Replace("$","\\$");
		rv.Replace("_","\\_");
		rv.Replace("#","\\#");
		return rv;
	}

	void BuildLatexEntry(mystr& entry, const LogData& data)
	{
		mystr date, dataentry, type;
		GetDate(date, data.date);
		entry = mystr(data.ver_num) + "." + mystr(data.rel_num) + "." + mystr(data.log_num) + " & ";
		entry += date;
		//entry += ", " + data.time
		entry += " & ";
		entry += AdaptStringToLatex(TransformUser(data.user)) + " & ";
		type = mystr((data.public_internal ? "p-" : "i-")) + data.type;
		entry += AdaptStringToLatex(TransformType(type));
		
		dataentry = AdaptStringToLatex(data.entry);
		//dataentry.Replace("&","\\&");
		entry += " & " + dataentry;

		//entry.Replace("\\","\\\\");
		//entry.Replace("$","\\$");
		//entry.Replace("_","\\_");
		//entry.Replace("#","\\#");
		entry.Replace("\n\n","@@@---@@@---@@@");
		entry.Replace("\n","\n\n");
		entry.Replace("@@@---@@@---@@@","\n\n\\vspace*{0.4cm}");
		entry += " \\\\ \\hline";
	}

	int WriteLatexTable(TArrayDynamic<LogData>* logdata)
	{
		int rv = 0;
		int opt = ios::out;
		//int opt = ios::out;
		ofstream os(latex_table_filename->c_str(), opt);
		mystr entry;
		
		if (os.is_open())
		{
			os << "\\begin{center}" << endl << "\\footnotesize" << endl;
			os << "  \\begin{longtable}{ | p{1.1cm} | p{1.7cm} | p{1.7cm} | p{1.4cm} | p{10.0cm} |} \\hline" << endl;
			os << "    \\bf log-\\# & \\bf date & \\bf user & \\bf type & \\bf description \\\\ \\hline" << endl;
			os << "    \\hline" << endl;;
			for (int i = 1; i <= logdata->Length(); ++i)
			{
				BuildLatexEntry(entry, (*logdata)(i));
				//int lognum = (*logdata)(i).log_num;
				os << "    " << entry << endl;
			}
			os <<	"  \\end{longtable}" << endl;
			os << "\\end{center}" << endl;
			os.close();
		}
		else
		{
			MessageBox::Show("The latex-file could not created!","Problem",MessageBoxButtons::OK,MessageBoxIcon::Error);
			rv = 1;
		}
		return rv;
	}

	int TransformLogToLatex()
	{
		int rv = 0;
		//int num_logs = DetermineNumOfLogs();
		//TArrayDynamic<LogData>* logdata = new TArrayDynamic<LogData>(num_logs);
		//rv = ReadLogs(logdata);
		rv += WriteLatexTable(logData);
		if (!rv)
		{
		}
		else
		{
			MessageBox::Show("Problem while reading logs!","Problem",MessageBoxButtons::OK,MessageBoxIcon::Error);
		}
		return rv;
	}

	int WriteHTMLFile(TArrayDynamic<LogData>* logdata, int public_internal)
	{
		int rv = 0;
		int opt = ios::out;
		//int opt = ios::out;
		ofstream *os;
		if (public_internal)
			os = new ofstream(html_public_filename->c_str(), opt);
		else
			os = new ofstream(html_filename->c_str(), opt);
		mystr entry;
		
		if (os->is_open())
		{
			(*os) << "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 3.2 Final//EN\">" << endl;
			(*os) << "<html><head><title>HOTINT Logs</title></head>" << endl;
			(*os) << "<body>" << endl;
			//(*os) << "<h3>HOTINT Logs</h3>" << endl;
			//(*os) << "<br><h4>Created by using HOTINT Logger</a></h4><p><table border=\"1\" cellpadding=\"5\"><tr bgcolor=\"E0E0E0\">" << endl;
			if (!public_internal)
				(*os) << "<br><h4>HOTINT Logs</a></h4><p><table border=\"1\" cellpadding=\"5\"><tr bgcolor=\"E0E0E0\">" << endl;
			else
				(*os) << "<br><h4>HOTINT Logs ignoring minor changes</a></h4><p><table border=\"1\" cellpadding=\"5\"><tr bgcolor=\"E0E0E0\">" << endl;
			(*os) << "<th>Log-Number" << endl;
			(*os) << "<th>Date" << endl;
			(*os) << "<th>User" << endl;
			(*os) << "<th>Flag" << endl;
			(*os) << "<th>Log-Message" << endl;
			for (int i = logdata->Length(); i >= 1; --i)
			{
				if ((!(*logdata)(i).public_internal && public_internal == 0) || (*logdata)(i).public_internal) 
				{	
					BuildHTMLEntry(entry, (*logdata)(i));
					//int lognum = (*logdata)(i).log_num;
					(*os) << entry << endl;
				}
			}
			(*os) <<	"</table>" << endl;
			(*os) << "<p><table border=\"0\" cellpadding=\"3\">" << endl;
			(*os) << "<th>type<tr><td><td>c<td>change" << endl;
			(*os) << "<tr><td><td>e<td>extension" << endl;
			(*os) << "<tr><td><td>f<td>new feature" << endl;
			(*os) << "<tr><td><td>b<td>bug fix" << endl;
			(*os) << "<tr><td><td>n<td>new" << endl;
			(*os) <<	"</table>" << endl;
			(*os) << "</body></html>" << endl;
			os->close();
		}
		else
		{
			MessageBox::Show("The HTML-table could not created! The html-file-path does not exist!","Problem",MessageBoxButtons::OK,MessageBoxIcon::Error);
			rv = 1;
		}
		delete os;
		return rv;
	}

	int WriteHTMLFile(TArrayDynamic<BugData>* bugdata, int public_internal)
	{
		int rv = 0;
		int opt = ios::out;
		//int opt = ios::out;
		ofstream *os;
		if (public_internal)
			os = new ofstream(html_bug_public_filename->c_str(), opt);
		else
			os = new ofstream(html_bug_filename->c_str(), opt);
		mystr entry;
		
		if (os->is_open())
		{
			(*os) << "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 3.2 Final//EN\">" << endl;
			(*os) << "<html><head><title>HOTINT Bugs</title></head>" << endl;
			(*os) << "<body>" << endl;
			//(*os) << "<h3>HOTINT Logs</h3>" << endl;
			//(*os) << "<br><h4>Created by using HOTINT Logger</a></h4><p><table border=\"1\" cellpadding=\"5\"><tr bgcolor=\"E0E0E0\">" << endl;
			if (!public_internal)
				(*os) << "<br><h4>HOTINT Bugs</a></h4><p><table border=\"1\" cellpadding=\"5\"><tr bgcolor=\"E0E0E0\">" << endl;
			else
				(*os) << "<br><h4>HOTINT Bugs ignoring minor changes</a></h4><p><table border=\"1\" cellpadding=\"5\"><tr bgcolor=\"E0E0E0\">" << endl;
			(*os) << "<th>Bug" << endl;
			(*os) << "<th>Date" << endl;
			(*os) << "<th>User" << endl;
			(*os) << "<th>Solved" << endl;
			(*os) << "<th>Bug-Description" << endl;
			(*os) << "<th>Bug-Solution" << endl;
			for (int i = bugdata->Length(); i >= 1; --i)
			{
				if ((!(*bugdata)(i).public_internal && public_internal == 0) || (*bugdata)(i).public_internal) 
				{	
					BuildHTMLEntry(entry, (*bugdata)(i));
					//int lognum = (*logdata)(i).log_num;
					(*os) << entry << endl;
				}
			}
			(*os) <<	"</table>" << endl;
			//(*os) << "<p><table border=\"0\" cellpadding=\"3\">" << endl;
			//(*os) << "<th>type<tr><td><td>c<td>change" << endl;
			//(*os) << "<tr><td><td>e<td>extension" << endl;
			//(*os) << "<tr><td><td>f<td>new feature" << endl;
			//(*os) << "<tr><td><td>b<td>bug fix" << endl;
			//(*os) << "<tr><td><td>n<td>new" << endl;
			//(*os) <<	"</table>" << endl;
			(*os) << "</body></html>" << endl;
			os->close();
		}
		else
		{
			MessageBox::Show("The HTML-table could not created! The html-file-path does not exist!","Problem",MessageBoxButtons::OK,MessageBoxIcon::Error);
			rv = 1;
		}
		delete os;
		return rv;
	}

	int TransformToHTML()
	{
		int rv = 0;
	  //int num_logs = DetermineNumOfLogs();
		//TArrayDynamic<LogData>* logdata = new TArrayDynamic<LogData>(num_logs);
		//rv = ReadLogs(logdata);
		if (isLog)
		{
			rv += WriteHTMLFile(logData, 0);
			rv += WriteHTMLFile(logData, 1);
		}
		else if (isBug)
		{
			rv += WriteHTMLFile(bugData, 0);
			rv += WriteHTMLFile(bugData, 1);
		}		
		if (!rv)
		{
			MessageBox::Show("Transformation to HTML was successful!","Transform to HTML",MessageBoxButtons::OK,MessageBoxIcon::Information );
		}
		else
		{
			MessageBox::Show("Problem in creating HTML table!","Problem",MessageBoxButtons::OK,MessageBoxIcon::Error);
		}
		return rv;
	}

	mystr TransformUser(const mystr& user)
	{
		//mystr rv(user);
		System::String^ rv2 = gcnew System::String(user.c_str());
		//System::String^ rv = rv2->ToLower();
		mystr rv;
		ConvertStringText(rv, rv2->ToLower());
		//System::String^ *str = gcnew System::String(rv.c_str());
		rv.Replace("schoergenhumer","sch");
		rv.Replace("dorninger","dor");
		rv.Replace("karer","kar");
		rv.Replace("gerstmayr","ger");
		rv.Replace("nachbagauer","nac");
		rv.Replace("reischl","rei");
		rv.Replace("ludwig","lud");
		rv.Replace("eder","ede");
		rv.Replace("gruber","gru");
		rv.Replace("vetyukov","vet");
		rv.Replace("pechstein","pec");
		rv.Replace("humer","hum");
		rv.Replace("stangl","sta");
		rv.Replace("saxinger","sax");
		rv.Replace("weitzhofer","wei");
		rv.Replace("aigner","aig");
		return rv;
	}

	mystr TransformType(const mystr& type)
	{
		mystr rv(type);
		rv.Replace("new features","f");
		rv.Replace("new feature","f");
		rv.Replace("critical change","cc");
		rv.Replace("change","c");
		rv.Replace("bugs fixed","b");
		rv.Replace("bug fixed","b");
		rv.Replace("bug fix","b");
		rv.Replace("bugfix","b");
		rv.Replace("extension","e");
		rv.Replace("new element","ne");
		rv.Replace("new","n");
		return rv;
	}

	mystr TransformHTML(const mystr& entry)
	{
		mystr rv(entry);
		rv.Replace("&","&amp;");
		rv.Replace("<","&lt;");
		rv.Replace(">","&gt;");
		rv.Replace("\"","&quot;");
		rv.Replace("´´","&quot;");
		rv.Replace(mystr('“'),"&quot;");
		return rv;
	}

	void BuildHTMLEntry(mystr& entry, const LogData& data)
	{
		mystr date;
		GetDate(date, data.date);
		/*entry = mystr(data.log_num) + " & ";
		entry += date + ", " + data.time + " & ";
		entry += data.user + " & " + data.type;
		entry += " & " + data.entry;

		entry.Replace("\\","\\\\");
		entry.Replace("$","\\$");
		entry.Replace("_","\\_");
		entry.Replace("\n\n","@@@---@@@---@@@");
		entry.Replace("\n","\n\n");
		entry.Replace("@@@---@@@---@@@","\n\n\\vspace*{0.4cm}");
		entry += " \\\\ \\hline";*/
		entry = "<tr><td bgcolor=#FFFFFF nowrap>" + mystr(data.ver_num) + "." + mystr(data.rel_num) + "." + mystr(data.log_num);
		entry += "<td bgcolor=#FFFFFF nowrap>" + date;
		entry += "<td bgcolor=#FFFFFF >" + TransformHTML(TransformUser(data.user));
		entry += "<td bgcolor=#FFFFFF >" + TransformHTML(TransformType(data.type));
		entry += "<td bgcolor=#FFFFFF >" + TransformHTML(data.entry);
		entry.Replace("\n","\n<br>");
	}

	void BuildHTMLEntry(mystr& entry, const BugData& data)
	{
		mystr date;
		GetDate(date, data.date);
		entry = "<tr><td bgcolor=#FFFFFF nowrap>" + mystr(data.ver_num) + "." + mystr(data.rel_num) + "." + mystr(data.log_num);
		entry += "<td bgcolor=#FFFFFF nowrap>" + date;
		entry += "<td bgcolor=#FFFFFF >" + TransformHTML(TransformUser(data.user));
		//entry += "<td bgcolor=#FFFFFF >" + TransformHTML(TransformType(data.type));
		entry += "<td bgcolor=#FFFFFF >" + TransformHTML((data.solved ? "yes" : "no"));
		entry += "<td bgcolor=#FFFFFF >" + TransformHTML(data.entry);
		entry += "<td bgcolor=#FFFFFF >" + TransformHTML(data.solution);
		entry.Replace("\n","\n<br>");
	}

	//in groupbox gBoxType the logs radio button was chosen -> logs will be displayed
	int ActivateLogs()
	{
		int rv = 0;
		rv += InitLogs();	
		rv += AdaptToLogFields();
		return rv;
	}

	//in groupbox gBoxType the bugs radio button was chosen -> bugs will be displayed
	int ActivateBugs()
	{
		int rv = 0;
		rv += InitBugs();
		rv += AdaptToBugFields();
		//bugChanged = false;
		BugChanged(false);
		return rv;
	}

	void BugChanged(bool status)
	{
		EnableBNav(!status);
		this->bNavAddNewItem->Enabled = !status;
		this->btnNew->Enabled = !status;
		this->btnTransformHTML->Enabled = !status;
		this->btnSave->Enabled = status;
		this->btnRefresh->Enabled = !status;
		this->btnCancel->Enabled = status;
		bugChanged = status;
	}

	int AdaptToBugFields()
	{
		int rv = 0;
		this->cbxBugSolved->Visible = true;
		this->lblNumEntry->Text = L"Bug-Number";
		this->lblLogDate->Text = L"Bug-Date/Time";
		this->lblLogType->Visible = false;
		this->cbLogType->Visible = false;
		this->lblBugSolution->Visible = true;
		this->rtbBugSolution->Visible = true;
		//this->cbxPublicInternal->Visible = true;
		this->rtbLogDescription->Size = System::Drawing::Size(643, 146);
			//System::Windows::Forms::Control::Size()
		return rv;
	}

	int AdaptToLogFields()
	{
		int rv = 0;
		this->cbxBugSolved->Visible = false;
		this->lblNumEntry->Text = L"Log-Number";
		this->lblLogDate->Text = L"Log-Date/Time";
		this->lblLogType->Visible = true;
		this->cbLogType->Visible = true;
		this->lblBugSolution->Visible = false;
		this->rtbBugSolution->Visible = false;
		//this->cbxPublicInternal->Visible = true;
		this->rtbLogDescription->Size = System::Drawing::Size(643, 224);
		return rv;
	}

	int InitBugs()
	{
		int rv = 0;

		int num_bugs = DetermineNumOfBugs();
		if (num_bugs > 0)
		{
			if (bugData)
				delete bugData;
			bugData = new TArrayDynamic<BugData>(num_bugs);
			rv += ReadBugs(bugData);
			//rv += LoadLastLog();
			rv += InitializeLogDataNavigation();
		}
			
		return rv;
	}


	System::Void btnTransform_Click(System::Object^  sender, System::EventArgs^  e) 
	{
		int rv = TransformStep1();
		rv += TransformStep2(); //add internal/public flag
		rv += InitializeLogDataNavigation();
		if (rv)
			MessageBox::Show("Transformation not successful!","Problem",MessageBoxButtons::OK,MessageBoxIcon::Error);
		else
			MessageBox::Show("Transformation successful!","Information",MessageBoxButtons::OK,MessageBoxIcon::Information);
	}
	System::Void btnTransformLatex_Click(System::Object^  sender, System::EventArgs^  e) 
	{
		int rv = TransformLogToLatex();
	}
  System::Void btnTransformHTML_Click(System::Object^  sender, System::EventArgs^  e) {
		int rv = TransformToHTML();	
	}
	System::Void bNavAddNewItem_Click(System::Object^  sender, System::EventArgs^  e) {
		btnNew_Click(sender, e);
	}
	System::Void bNavMoveLastItem_Click(System::Object^  sender, System::EventArgs^  e) {
		//ShowLog(logData->Last());
		if (isLog)
			this->bNavPositionItem->Text = System::Convert::ToString(logData->Length());
		else if (isBug)
			this->bNavPositionItem->Text = System::Convert::ToString(bugData->Length());
		//this->bNavMoveLastItem->Enabled = false; 
		//this->bNavMoveFirstItem->Enabled = true; 
		//this->bNavPositionItem->Enabled = true;
	}
	System::Void bNavMoveNextItem_Click(System::Object^  sender, System::EventArgs^  e) {
		int num_entry = System::Convert::ToInt32(this->bNavPositionItem->Text);
		if (isLog)
			num_entry = num_entry%logData->Length() + 1;
		else if (isBug)
			num_entry = num_entry%bugData->Length() + 1;
		//ShowLog((*logData)(log));
		this->bNavPositionItem->Text = System::Convert::ToString(num_entry);
	}
	System::Void bNavMovePreviousItem_Click(System::Object^  sender, System::EventArgs^  e) {
		int num_entry = System::Convert::ToInt32(this->bNavPositionItem->Text);
		if (isLog)
			num_entry = (num_entry + logData->Length() - 2)%logData->Length() + 1;
		else if (isBug)
			num_entry = (num_entry + bugData->Length() - 2)%bugData->Length() + 1;
		this->bNavPositionItem->Text = System::Convert::ToString(num_entry);
	}
	System::Void bNavMoveFirstItem_Click(System::Object^  sender, System::EventArgs^  e) {
		this->bNavPositionItem->Text = System::Convert::ToString(1);
	}

	System::Void bNavPositionItem_TextChanged(System::Object^  sender, System::EventArgs^  e) 
	{
		if (isNew) return;
		int num_entry;
		if (isBug && bugChanged) 
		{
			try
			{
				num_entry = System::Convert::ToInt32(this->bNavPositionItem->Text);
				if (num_entry != curBug)
					this->bNavPositionItem->Text = System::Convert::ToString(curBug); 
			}
			catch(System::Exception^)		
			{
				this->bNavPositionItem->Text = System::Convert::ToString(curBug);
			}
			return;
		}
		try
		{
			num_entry = System::Convert::ToInt32(this->bNavPositionItem->Text);
		}
		catch(System::Exception^)		
		{
			if ((!logData && isLog) || (!bugData && isBug))
			{
				this->bNavPositionItem->Text = System::Convert::ToString(0);
				return;
			}	
			else if (isLog)
			{	num_entry = logData->Length(); }
			else if (isBug)
			{ num_entry = bugData->Length(); }
		}
		if ((!logData && isLog) || (!bugData && isBug))
			return;
		int len = (isLog ? logData->Length() : (isBug ? bugData->Length() : 0));
		num_entry = (num_entry-1)%len + 1;
		if (num_entry <= 0)
			num_entry = len;
		if (isLog)
			ShowLog((*logData)(num_entry));
		else if (isBug)
		{
			ShowBug((*bugData)(num_entry));
			curBug = num_entry;
		}
		this->bNavPositionItem->Text = System::Convert::ToString(num_entry);
		this->bNavMoveLastItem->Enabled = (num_entry == len) ? false : true; 
		this->bNavMoveFirstItem->Enabled = (num_entry == 1) ? false : true; 
		this->bNavPositionItem->Enabled = true;
	}
	
	System::Void btnCancel_Click(System::Object^  sender, System::EventArgs^  e) {
		if (isLog && !isNew) return;
		EnableButtons(true);
		EnableFields(false);
		if (isLog) // btnCancel should only be available if isNew & isLog
		{
			InitializeLogDataNavigation();
			LoadLastLog();
		}
		else if (isBug)
		{
			if (isNew)
			{
				InitializeLogDataNavigation();
				ShowBug(bugData->Last());
				curBug = bugData->Length();
			}
			else
			{
				EnableBNav(true);
				curBug = System::Convert::ToInt32(this->bNavPositionItem->Text);
				ShowBug((*bugData)(curBug));
				if (curBug == 1)
					this->bNavMoveFirstItem->Enabled = false;
				if (curBug == bugData->Length())
					this->bNavMoveLastItem->Enabled = false;
			}
		}
		isNew = false;
	}

	System::Void btnRefresh_Click(System::Object^  sender, System::EventArgs^  e) {
	  Destroy();
		int identType = (isLog ? 1 : (isBug ? 2 : -1));
		int rv = Initialize(identType);
	}
	System::Void rbLog_CheckedChanged(System::Object^  sender, System::EventArgs^  e) 
	{
		int rv = 0;
		isLog = rbLog->Checked;
		isBug = rbBug->Checked;
		if (rbLog->Checked)
		{	rv += ActivateLogs(); }
		else if (rbBug->Checked)
		{	rv += ActivateBugs(); }
		else if (rbToDo->Checked)
		{
		}
		else
		{
		}
	}
	System::Void cbxBugSolved_CheckedChanged(System::Object^  sender, System::EventArgs^  e) 
	{
		BugChanged(true);
	}
  System::Void rtbBugSolution_TextChanged(System::Object^  sender, System::EventArgs^  e) 
	{
		BugChanged(true);
	}



};


	/*struct Date
	{	
		int day; 
		int month; 
		int year; 
	
		Date& operator= (Date& date)
		{
			day = date.day;
			month = date.month;
			year = date.year;
			return date;
		}
	};*/


//	class BaseData
//	{
//	public:
//		int ver_num;
//		int rel_num;
//		int log_num;
//		Date date;
//		mystr time;
//		mystr user;									
//		mystr entry;								//the log entry / description
//			
//		BaseData() : ver_num(0), rel_num(0), log_num(0), date(), time(), user(), entry() 
//		{}
//
//		BaseData(BaseData& data)
//		{
//			CopyFrom(data);
//		}
//
//		void CopyFrom(BaseData& data)
//		{
//			ver_num = data.ver_num;
//			rel_num = data.rel_num;
//			log_num = data.log_num;
//			date = data.date;
//			time = data.time;
//			user = data.user;
//			entry = data.entry;
//		}
//	};
//	
//
//
//	class LogData : public BaseData
//	{
//	public:
//		mystr type;									//the log type ... new feature (f), change (c), extension (e)
//		bool public_internal;				//defines whether a log is public or internal
//	
//		LogData() : BaseData()
//		{ }
//
//		LogData(LogData& data)
//		{
//			CopyFrom(data);
//		}
//
//		void CopyFrom(LogData& data)
//		{
//			BaseData::CopyFrom(data);
//			type = data.type;
//			public_internal = data.public_internal;
//		}
//	};
//
//	////struct defining the data of the log
//	//struct LogData 
//	//{
//	//	int ver_num;
//	//	int rel_num;
//	//	int log_num;
//	//	Date date;
//	//	mystr time;
//	//	mystr user;									
//	//	mystr type;									//the log type ... new feature (f), change (c), extension (e)
//	//	mystr entry;								//the log entry / description
//	//	bool public_internal;				//defines whether a log is public or internal
//	//};
//
//	class BugData: public BaseData
//	{
//	public:
//		mystr solution;							//possibly the solution to the bug....
//		bool solved;								//bug solved??
//	
//		BugData() : BaseData()
//		{ }
//
//		BugData(BugData& data)
//		{
//			CopyFrom(data);
//		}
//
//		void CopyFrom(BugData& data)
//		{
//			BaseData::CopyFrom(data);
//			solution = data.solution;
//			solved = data.solved;
//		}
//	};
//
//	//struct BugData 
//	//{
//	//	int ver_num;
//	//	int rel_num;
//	//	int bug_num;
//	//	Date date;
//	//	mystr time;
//	//	mystr user;									
//	//	mystr entry;								//the log entry / description
//	//	mystr solution;							//possibly the solution to the bug....
//	//	bool solved;								//bug solved??	
//	//};
//
//
//	//Sorts an array a(1..a.Length()) into ascending numerical order by Shell’s method (diminishing increment
////sort). a is replaced on output by its sorted rearrangement.
//
//template <typename T, typename T2>
//void Quicksort(TArray<T>& a, TArray<T2>& b)
//{
//	int n = a.Length();
//	int i,j,inc;
//	T v;
//	T2 vb;
//	inc=1; //Determine the starting increment.
//	do 
//	{
//		inc *= 3;
//		inc++;
//	} while (inc <= n);
//	do 
//	{ //Loop over the partial sorts.
//		inc /= 3;
//		for (i=inc;i<n;i++) 
//		{ //Outer loop of straight insertion.
//			v = a.Get0(i);
//			vb = b.Get0(i);
//			j=i;
//			while (a.Get0(j-inc) > v) 
//			{ //Inner loop of straight insertion.
//				a.Elem0(j) = a.Get0(j-inc);
//				b.Elem0(j) = b.Get0(j-inc);
//				j -= inc;
//				if (j < inc) break;
//			}
//			a.Elem0(j)=v;
//			b.Elem0(j)=vb;
//		}
//	} while (inc > 1);
//
//}
//
//
//

}

