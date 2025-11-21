/**
 * Project List Page - Complete project management interface
 */
import React, { useEffect, useState } from 'react'
import {
  Card,
  Button,
  Table,
  Tag,
  Space,
  Modal,
  Dropdown,
  Statistic,
  Row,
  Col,
  Input,
  Select,
  message as antMessage,
} from 'antd'
import {
  PlusOutlined,
  FolderOutlined,
  EditOutlined,
  DeleteOutlined,
  MoreOutlined,
  InboxOutlined,
  RestOutlined,
  CheckCircleOutlined,
  ClockCircleOutlined,
  FolderOpenOutlined,
} from '@ant-design/icons'
import type { ColumnsType } from 'antd/es/table'
import { useProjectStore } from '../../store/projectStore'
import { ProjectFormModal } from './components/ProjectFormModal'
import type { Project } from '../../types/project'
import dayjs from 'dayjs'
import relativeTime from 'dayjs/plugin/relativeTime'

dayjs.extend(relativeTime)

const { Search } = Input
const { Option } = Select

export const ProjectList: React.FC = () => {
  const {
    projects,
    stats,
    loading,
    fetchProjects,
    fetchStats,
    deleteProject,
    archiveProject,
    restoreProject,
  } = useProjectStore()

  const [modalOpen, setModalOpen] = useState(false)
  const [editingProject, setEditingProject] = useState<Project | null>(null)
  const [statusFilter, setStatusFilter] = useState<string>('all')
  const [searchText, setSearchText] = useState('')

  useEffect(() => {
    fetchProjects()
    fetchStats()
  }, [])

  // Filter projects based on status and search
  const filteredProjects = projects.filter((project) => {
    const matchesStatus = statusFilter === 'all' || project.status === statusFilter
    const matchesSearch =
      searchText === '' ||
      project.name.toLowerCase().includes(searchText.toLowerCase()) ||
      project.description?.toLowerCase().includes(searchText.toLowerCase())
    return matchesStatus && matchesSearch
  })

  const handleCreate = () => {
    setEditingProject(null)
    setModalOpen(true)
  }

  const handleEdit = (project: Project) => {
    setEditingProject(project)
    setModalOpen(true)
  }

  const handleDelete = (project: Project) => {
    Modal.confirm({
      title: 'Delete Project',
      content: `Are you sure you want to delete "${project.name}"? This will also delete all associated samples, files, and tasks.`,
      okText: 'Delete',
      okType: 'danger',
      onOk: async () => {
        await deleteProject(project.id)
      },
    })
  }

  const handleArchive = async (project: Project) => {
    await archiveProject(project.id)
  }

  const handleRestore = async (project: Project) => {
    await restoreProject(project.id)
  }

  const columns: ColumnsType<Project> = [
    {
      title: 'Project Name',
      dataIndex: 'name',
      key: 'name',
      render: (name, record) => (
        <Space>
          <FolderOutlined style={{ fontSize: 18, color: '#1890ff' }} />
          <span style={{ fontWeight: 500 }}>{name}</span>
          {record.status === 'archived' && <Tag color="default">Archived</Tag>}
        </Space>
      ),
    },
    {
      title: 'Description',
      dataIndex: 'description',
      key: 'description',
      ellipsis: true,
    },
    {
      title: 'Status',
      dataIndex: 'status',
      key: 'status',
      width: 120,
      render: (status) => {
        const statusConfig: Record<string, { color: string; icon: React.ReactNode }> = {
          active: { color: 'success', icon: <CheckCircleOutlined /> },
          archived: { color: 'default', icon: <InboxOutlined /> },
          completed: { color: 'blue', icon: <FolderOpenOutlined /> },
        }
        const config = statusConfig[status] || statusConfig.active
        return (
          <Tag color={config.color} icon={config.icon}>
            {status.toUpperCase()}
          </Tag>
        )
      },
    },
    {
      title: 'Samples',
      dataIndex: 'sample_count',
      key: 'sample_count',
      width: 100,
      render: (count) => count || 0,
    },
    {
      title: 'Tasks',
      dataIndex: 'task_count',
      key: 'task_count',
      width: 100,
      render: (count) => count || 0,
    },
    {
      title: 'Created',
      dataIndex: 'created_at',
      key: 'created_at',
      width: 150,
      render: (date) => dayjs(date).fromNow(),
    },
    {
      title: 'Actions',
      key: 'actions',
      width: 100,
      render: (_, record) => (
        <Dropdown
          menu={{
            items: [
              {
                key: 'edit',
                label: 'Edit',
                icon: <EditOutlined />,
                onClick: () => handleEdit(record),
              },
              {
                key: 'archive',
                label: record.status === 'archived' ? 'Restore' : 'Archive',
                icon: record.status === 'archived' ? <RestOutlined /> : <InboxOutlined />,
                onClick: () =>
                  record.status === 'archived'
                    ? handleRestore(record)
                    : handleArchive(record),
              },
              {
                type: 'divider',
              },
              {
                key: 'delete',
                label: 'Delete',
                icon: <DeleteOutlined />,
                danger: true,
                onClick: () => handleDelete(record),
              },
            ],
          }}
        >
          <Button icon={<MoreOutlined />} />
        </Dropdown>
      ),
    },
  ]

  return (
    <div>
      {/* Statistics Cards */}
      <Row gutter={16} style={{ marginBottom: 24 }}>
        <Col xs={24} sm={12} lg={6}>
          <Card>
            <Statistic
              title="Total Projects"
              value={stats?.total_projects || 0}
              prefix={<FolderOutlined />}
              valueStyle={{ color: '#1890ff' }}
            />
          </Card>
        </Col>
        <Col xs={24} sm={12} lg={6}>
          <Card>
            <Statistic
              title="Active Projects"
              value={stats?.active_projects || 0}
              prefix={<CheckCircleOutlined />}
              valueStyle={{ color: '#52c41a' }}
            />
          </Card>
        </Col>
        <Col xs={24} sm={12} lg={6}>
          <Card>
            <Statistic
              title="Total Tasks"
              value={stats?.total_tasks || 0}
              prefix={<ClockCircleOutlined />}
              valueStyle={{ color: '#722ed1' }}
            />
          </Card>
        </Col>
        <Col xs={24} sm={12} lg={6}>
          <Card>
            <Statistic
              title="Active Tasks"
              value={stats?.active_tasks || 0}
              prefix={<ClockCircleOutlined />}
              valueStyle={{ color: '#fa8c16' }}
            />
          </Card>
        </Col>
      </Row>

      {/* Filters and Actions */}
      <Card style={{ marginBottom: 16 }}>
        <Space style={{ width: '100%', justifyContent: 'space-between' }}>
          <Space>
            <Search
              placeholder="Search projects..."
              style={{ width: 300 }}
              value={searchText}
              onChange={(e) => setSearchText(e.target.value)}
              allowClear
            />
            <Select
              value={statusFilter}
              onChange={setStatusFilter}
              style={{ width: 150 }}
            >
              <Option value="all">All Status</Option>
              <Option value="active">Active</Option>
              <Option value="archived">Archived</Option>
              <Option value="completed">Completed</Option>
            </Select>
          </Space>
          <Button type="primary" icon={<PlusOutlined />} onClick={handleCreate}>
            New Project
          </Button>
        </Space>
      </Card>

      {/* Projects Table */}
      <Card>
        <Table
          columns={columns}
          dataSource={filteredProjects}
          rowKey="id"
          loading={loading}
          pagination={{
            pageSize: 10,
            showSizeChanger: true,
            showTotal: (total) => `Total ${total} projects`,
          }}
        />
      </Card>

      {/* Create/Edit Modal */}
      <ProjectFormModal
        open={modalOpen}
        project={editingProject}
        onClose={() => {
          setModalOpen(false)
          setEditingProject(null)
        }}
        onSuccess={() => {
          setModalOpen(false)
          setEditingProject(null)
          fetchProjects()
          fetchStats()
        }}
      />
    </div>
  )
}

export default ProjectList
